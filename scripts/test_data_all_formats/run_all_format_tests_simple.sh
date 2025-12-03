#!/bin/bash
################################################################################
# scRNA-DataHub 全格式测试脚本 (简化版)
# 不使用 set -e，手动检查每个命令的返回值
################################################################################

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

PASSED=0
FAILED=0
SKIPPED=0

print_header() {
    echo -e "\n${BLUE}========================================================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================================================${NC}\n"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_info() {
    echo -e "${BLUE}ℹ $1${NC}"
}

test_format() {
    local format_name="$1"
    local input_path="$2"
    local output_file="$3"
    local extra_args="$4"
    
    echo -e "\n${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${YELLOW}测试 #$((PASSED + FAILED + SKIPPED + 1)): $format_name${NC}"
    echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    
    # 检查输入文件是否存在
    if [ ! -e "$input_path" ]; then
        print_warning "输入文件不存在，跳过"
        echo -e "  ${YELLOW}文件路径: $input_path${NC}"
        ((SKIPPED++))
        echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
        return 0
    fi
    
    echo -e "  📂 ${BLUE}输入文件: $input_path${NC}"
    echo -e "  💾 ${BLUE}输出文件: $output_file${NC}"
    echo ""
    
    # 运行测试，捕获输出
    local temp_output=$(mktemp)
    python "$READER_SCRIPT" "$input_path" -o "$output_file" $extra_args 2>&1 | tee "$temp_output"
    local test_result=${PIPESTATUS[0]}
    
    echo ""
    
    if [ $test_result -eq 0 ]; then
        # 验证输出文件
        if [ -f "$output_file" ]; then
            local file_size=$(du -h "$output_file" | cut -f1)
            
            # 读取H5AD文件统计信息
            local stats=$(python -c "
import scanpy as sc
import sys
try:
    adata = sc.read_h5ad('$output_file')
    print(f'CELLS:{adata.n_obs}')
    print(f'GENES:{adata.n_vars}')
    print(f'SPARSE:{\"是\" if hasattr(adata.X, \"toarray\") else \"否\"}')
    obs_cols = list(adata.obs.columns) if len(adata.obs.columns) > 0 else ['无']
    var_cols = list(adata.var.columns) if len(adata.var.columns) > 0 else ['无']
    print(f'OBS_COLS:{\" \".join(obs_cols[:5])}')
    print(f'VAR_COLS:{\" \".join(var_cols[:5])}')
except Exception as e:
    print(f'ERROR:{e}', file=sys.stderr)
    sys.exit(1)
" 2>&1)
            
            if [ $? -eq 0 ]; then
                # 解析统计信息
                local n_cells=$(echo "$stats" | grep "^CELLS:" | cut -d: -f2)
                local n_genes=$(echo "$stats" | grep "^GENES:" | cut -d: -f2)
                local is_sparse=$(echo "$stats" | grep "^SPARSE:" | cut -d: -f2)
                local obs_cols=$(echo "$stats" | grep "^OBS_COLS:" | cut -d: -f2-)
                local var_cols=$(echo "$stats" | grep "^VAR_COLS:" | cut -d: -f2-)
                
                print_success "✅ $format_name 测试通过"
                echo ""
                echo -e "  ${GREEN}📊 数据统计:${NC}"
                echo -e "     • 细胞数: ${GREEN}$n_cells${NC}"
                echo -e "     • 基因数: ${GREEN}$n_genes${NC}"
                echo -e "     • 文件大小: ${GREEN}$file_size${NC}"
                echo -e "     • 稀疏矩阵: ${GREEN}$is_sparse${NC}"
                
                if [ "$obs_cols" != "无" ] && [ -n "$obs_cols" ]; then
                    echo -e "     • obs列: ${GREEN}$obs_cols${NC}"
                fi
                if [ "$var_cols" != "无" ] && [ -n "$var_cols" ]; then
                    echo -e "     • var列: ${GREEN}$var_cols${NC}"
                fi
                
                ((PASSED++))
                
                # 记录成功的测试
                echo "$format_name|$n_cells|$n_genes|$file_size|$output_file" >> "$TEST_LOG"
            else
                print_success "$format_name 文件生成成功 (文件大小: $file_size)"
                print_warning "无法读取统计信息"
                ((PASSED++))
                echo "$format_name|N/A|N/A|$file_size|$output_file" >> "$TEST_LOG"
            fi
        else
            print_error "❌ $format_name 测试失败: 输出文件未生成"
            echo -e "  ${RED}预期输出: $output_file${NC}"
            ((FAILED++))
            echo "$format_name|FAILED|文件未生成|N/A|$output_file" >> "$FAIL_LOG"
        fi
    else
        print_error "❌ $format_name 测试失败 (退出码: $test_result)"
        echo -e "  ${RED}输入文件: $input_path${NC}"
        echo -e "  ${RED}输出文件: $output_file${NC}"
        
        # 提取错误信息
        local error_msg=$(grep -E "(Error|错误|Traceback)" "$temp_output" | tail -5)
        if [ -n "$error_msg" ]; then
            echo -e "  ${RED}错误信息:${NC}"
            echo "$error_msg" | while read line; do
                echo -e "    ${RED}$line${NC}"
            done
        fi
        
        ((FAILED++))
        echo "$format_name|FAILED|退出码$test_result|N/A|$input_path" >> "$FAIL_LOG"
    fi
    
    rm -f "$temp_output"
    echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    
    return 0  # 总是返回成功，继续下一个测试
}

################################################################################
# 主程序
################################################################################

print_header "scRNA-DataHub 全格式测试 (简化版)"

# 1. 获取测试数据目录路径
echo -e "${BLUE}请输入测试数据目录的绝对路径:${NC}"
echo -e "${YELLOW}示例: /storeData/ztron/wangrm/tools/scRNA-DataHub/scripts/test_data_all_formats${NC}"
read -p "路径: " TEST_DATA_DIR

# 验证路径
if [ ! -d "$TEST_DATA_DIR" ]; then
    print_error "目录不存在: $TEST_DATA_DIR"
    exit 1
fi

TEST_DATA_DIR=$(cd "$TEST_DATA_DIR" && pwd)
print_success "测试数据目录: $TEST_DATA_DIR"

# 2. 定位reader脚本
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
READER_SCRIPT="$PROJECT_ROOT/src/universal_reader.py"

if [ ! -f "$READER_SCRIPT" ]; then
    print_error "找不到reader脚本: $READER_SCRIPT"
    exit 1
fi
print_success "Reader脚本: $READER_SCRIPT"

# 3. 创建输出目录
OUTPUT_DIR="$TEST_DATA_DIR/test_outputs"
mkdir -p "$OUTPUT_DIR"
print_success "输出目录: $OUTPUT_DIR"

# 4. 清理旧的输出文件
print_info "清理旧的测试输出..."
rm -f "$OUTPUT_DIR"/*.h5ad 2>/dev/null || true
print_success "清理完成"

# 5. 创建日志文件
TEST_LOG="$OUTPUT_DIR/test_results.log"
FAIL_LOG="$OUTPUT_DIR/test_failures.log"
SUMMARY_FILE="$OUTPUT_DIR/test_summary.txt"

# 清空旧日志
> "$TEST_LOG"
> "$FAIL_LOG"

print_header "开始格式测试"
echo "测试时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "日志文件: $TEST_LOG"
echo ""

# 开始测试
test_format "10X MTX v3 (压缩)" "$TEST_DATA_DIR/filtered_feature_bc_matrix" "$OUTPUT_DIR/10x_mtx_v3.h5ad" "--sample-id 10x_mtx_v3"
test_format "10X MTX v2 (未压缩)" "$TEST_DATA_DIR/hg19_chr21" "$OUTPUT_DIR/10x_mtx_v2.h5ad" "--sample-id 10x_mtx_v2 --no-compressed"
test_format "10X H5 v3" "$TEST_DATA_DIR/filtered_feature_bc_matrix.h5" "$OUTPUT_DIR/10x_h5_v3.h5ad" "--sample-id 10x_h5_v3"
test_format "10X H5 v2" "$TEST_DATA_DIR/filtered_gene_bc_matrices_h5.h5" "$OUTPUT_DIR/10x_h5_v2.h5ad" "--sample-id 10x_h5_v2"
test_format "H5AD" "$TEST_DATA_DIR/10x_pbmc68k_reduced.h5ad" "$OUTPUT_DIR/h5ad.h5ad" "--sample-id h5ad"
test_format "Loom" "$TEST_DATA_DIR/test_data.loom" "$OUTPUT_DIR/loom.h5ad" "--sample-id loom"
test_format "Zarr" "$TEST_DATA_DIR/10x-10k-subset.zarr" "$OUTPUT_DIR/zarr.h5ad" "--sample-id zarr"
test_format "CSV" "$TEST_DATA_DIR/test_expression.csv" "$OUTPUT_DIR/csv.h5ad" "--sample-id csv --delimiter ,"
test_format "CSV (压缩)" "$TEST_DATA_DIR/test_expression.csv.gz" "$OUTPUT_DIR/csv_gz.h5ad" "--sample-id csv_gz --delimiter ,"
test_format "TSV" "$TEST_DATA_DIR/test_expression.tsv" "$OUTPUT_DIR/tsv.h5ad" "--sample-id tsv"
test_format "TSV (压缩)" "$TEST_DATA_DIR/test_expression.tsv.gz" "$OUTPUT_DIR/tsv_gz.h5ad" "--sample-id tsv_gz"
test_format "TXT" "$TEST_DATA_DIR/krumsiek11.txt" "$OUTPUT_DIR/txt.h5ad" "--sample-id txt"
test_format "Excel" "$TEST_DATA_DIR/test_expression.xlsx" "$OUTPUT_DIR/excel.h5ad" "--sample-id excel --sheet 0"
test_format "MTX单文件" "$TEST_DATA_DIR/test_matrix.mtx" "$OUTPUT_DIR/mtx.h5ad" "--sample-id mtx"
test_format "MTX单文件 (压缩)" "$TEST_DATA_DIR/test_matrix.mtx.gz" "$OUTPUT_DIR/mtx_gz.h5ad" "--sample-id mtx_gz"
test_format "UMI-tools" "$TEST_DATA_DIR/umi_tools_counts.tsv.gz" "$OUTPUT_DIR/umi_tools.h5ad" "--sample-id umi_tools"
test_format "自定义10X MTX" "$TEST_DATA_DIR/custom_10x_mtx" "$OUTPUT_DIR/custom_10x_mtx.h5ad" "--sample-id custom_10x_mtx"

################################################################################
# 测试结果汇总
################################################################################

print_header "测试结果汇总"

TOTAL=$((PASSED + FAILED + SKIPPED))
END_TIME=$(date '+%Y-%m-%d %H:%M:%S')

echo "测试完成时间: $END_TIME"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo -e "${BLUE}总体统计${NC}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  总测试数: $TOTAL"
echo -e "  ${GREEN}✅ 通过: $PASSED${NC}"
echo -e "  ${RED}❌ 失败: $FAILED${NC}"
echo -e "  ${YELLOW}⊘ 跳过: $SKIPPED${NC}"
echo -e "  通过率: $(awk "BEGIN {printf \"%.1f%%\", ($PASSED/$TOTAL)*100}")"
echo ""

# 显示成功的测试详情
if [ $PASSED -gt 0 ] && [ -f "$TEST_LOG" ]; then
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo -e "${GREEN}✅ 通过的测试 ($PASSED个)${NC}"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    printf "%-25s %10s %10s %10s\n" "格式" "细胞数" "基因数" "文件大小"
    echo "────────────────────────────────────────────────────────────────────"
    
    while IFS='|' read -r format cells genes size output; do
        if [ "$cells" != "N/A" ]; then
            printf "%-25s %10s %10s %10s\n" "$format" "$cells" "$genes" "$size"
        else
            printf "%-25s %10s %10s %10s\n" "$format" "-" "-" "$size"
        fi
    done < "$TEST_LOG"
    echo ""
fi

# 显示失败的测试详情
if [ $FAILED -gt 0 ] && [ -f "$FAIL_LOG" ]; then
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo -e "${RED}❌ 失败的测试 ($FAILED个)${NC}"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    
    local fail_count=1
    while IFS='|' read -r format status reason size path; do
        echo -e "${RED}$fail_count. $format${NC}"
        echo "   原因: $reason"
        echo "   文件: $path"
        echo ""
        ((fail_count++))
    done < "$FAIL_LOG"
fi

# 生成汇总报告文件
cat > "$SUMMARY_FILE" << EOF
========================================================================
scRNA-DataHub 测试报告
========================================================================

测试时间: $END_TIME
测试目录: $TEST_DATA_DIR

总体统计:
  - 总测试数: $TOTAL
  - 通过: $PASSED
  - 失败: $FAILED
  - 跳过: $SKIPPED
  - 通过率: $(awk "BEGIN {printf \"%.1f%%\", ($PASSED/$TOTAL)*100}")

详细结果:

通过的测试 ($PASSED个):
$(cat "$TEST_LOG" 2>/dev/null | awk -F'|' '{printf "  ✅ %-25s 细胞: %-8s 基因: %-8s 大小: %s\n", $1, $2, $3, $4}')

失败的测试 ($FAILED个):
$(cat "$FAIL_LOG" 2>/dev/null | awk -F'|' '{printf "  ❌ %-25s %s\n", $1, $3}')

输出文件位置: $OUTPUT_DIR
详细日志: $TEST_LOG
失败日志: $FAIL_LOG

========================================================================
EOF

print_success "测试报告已保存到: $SUMMARY_FILE"

# 显示输出文件列表
if [ $PASSED -gt 0 ]; then
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo -e "${BLUE}输出文件列表${NC}"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    ls -lh "$OUTPUT_DIR"/*.h5ad 2>/dev/null | awk '{printf "  %s %s %s\n", $9, $5, $6" "$7" "$8}' || echo "  (无输出文件)"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ $FAILED -eq 0 ]; then
    print_success "🎉 所有测试通过！"
    echo ""
    echo "下一步:"
    echo "  1. 查看详细报告: cat $SUMMARY_FILE"
    echo "  2. 验证输出文件: ls -lh $OUTPUT_DIR/*.h5ad"
    echo "  3. 读取数据示例: python -c 'import scanpy as sc; adata = sc.read_h5ad(\"$OUTPUT_DIR/10x_mtx_v3.h5ad\"); print(adata)'"
    exit 0
else
    print_error "⚠️  有 $FAILED 个测试失败"
    echo ""
    echo "故障排除:"
    echo "  1. 查看失败详情: cat $FAIL_LOG"
    echo "  2. 查看完整报告: cat $SUMMARY_FILE"
    echo "  3. 查看文档: cat TROUBLESHOOTING.md"
    exit 1
fi

