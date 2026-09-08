#!/bin/bash

# 设置脚本在遇到错误时立即退出
set -e

# --- 配置 ---
# 这个路径是脚本在 *容器内部* 看到的数据路径
BASE_DIR="/data"

# --- 主逻辑 ---
echo "=========================================="
echo "开始为所有 *.sorted.gff.gz 创建 Trix 文本索引..."
echo "每个索引将被创建在其对应的物种文件夹下。"
echo "=========================================="

# 检查基础目录是否存在
if [ ! -d "$BASE_DIR" ]; then
    echo "错误：容器内的目录 $BASE_DIR 不存在！请检查 docker -v 挂载参数。"
    exit 1
fi

# 使用 find 命令查找所有已经排序和压缩好的 GFF 文件
find "$BASE_DIR" -type f -name "*.sorted.gff.gz" | while IFS= read -r gff_gz_file; do

    # 获取包含该文件的目录，例如 /data/818.1357
    species_dir=$(dirname "$gff_gz_file")
    
    # 从目录路径中提取 genome_id
    genome_id=$(basename "$species_dir")
    
    echo "--- 正在为基因组: $genome_id 创建索引 ---"
    echo "输入文件: $gff_gz_file"
    
    # --- 执行带有新参数的 jbrowse text-index 命令 ---
    # --out 参数指向物种目录，jbrowse 会在此目录下创建 trix 子目录
    npx jbrowse text-index \
        --file "$gff_gz_file" \
        --out "$species_dir/" \
        --attributes "gene" \
        --fileId "${genome_id}_features" \
        --exclude "None" \
        --force
        
    echo "索引成功创建于: ${species_dir}/trix/"
    echo ""
done

echo "=========================================="
echo "所有 Trix 索引创建完毕！"
echo "=========================================="