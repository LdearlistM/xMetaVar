#!/bin/bash

# --- 配置 ---
# 包含所有待处理 fna.gz 文件的根目录
FASTA_DIR="/home/data/CYM/Webserver/xMetaVar/src/public/jbrowse_data/GENOMES_SV"

# --- 检查目录是否存在 ---
if [ ! -d "$FASTA_DIR" ]; then
    echo "错误：FASTA 目录 $FASTA_DIR 未找到！"
    exit 1
fi

# --- 主逻辑 ---
echo "开始处理 FASTA (.fna.gz) 文件..."

# 查找 FASTA_DIR 下所有的 .fna.gz 文件
find "$FASTA_DIR" -type f -name "*.fna.gz" | while IFS= read -r fasta_file; do
    echo "正在处理: $fasta_file"

    # 从文件路径中提取 genome ID (即其父目录的名称)
    genome_id=$(basename "$(dirname "$fasta_file")")
    
    # 定义临时输出文件名，这次我们直接输出未压缩的临时文件
    temp_output_file="${fasta_file}.tmp"

    # --- 核心处理逻辑 ---
    # 1. zcat: 解压文件并将内容输出到标准输出
    # 2. awk: 合并序列，只保留第一个 header
    # 3. sed: 修改保留下来的那个 header
    # 4. > temp_output_file: 将最终结果重定向到临时文件
    zcat "$fasta_file" | \
    awk '
        # 对以 ">" 开头的行 (header)
        /^>/ {
            # 如果是第一个 header (NR==1), 就打印它
            if (NR==1) {
                print $0
            }
            # 对于其他 header，直接跳过 (不打印)
            next
        }
        # 对于不是 header 的行 (序列行)
        {
            # 使用 printf 而不是 print，避免在行尾自动添加换行符
            # 这样所有的序列行就连成了一整行
            printf "%s", $0
        }
        # END 块在所有行处理完后执行
        END {
            # 在所有序列合并完后，打印一个换行符，格式更美观
            print ""
        }
    ' | \
    sed "1s/^>[^ ]\+/>$genome_id/" > "$temp_output_file"
    # sed "1s/^>[^ ]\+/>$genome_id/" 的解释:
    # 1s: 只对第一行进行操作
    # /^>[^ ]\+/: 匹配以 ">" 开头，后面跟着一个或多个非空格字符的模式 (即 >NC_015067)
    # />$genome_id/: 将匹配到的内容替换为 ">" 加上我们提取的 genome_id

    # 检查临时文件是否创建成功且非空
    if [ -s "$temp_output_file" ]; then
        # 将处理后的临时文件压缩并覆盖原始文件
        # gzip -c 会将压缩结果输出到标准输出，我们再重定向覆盖原文件
        gzip -c "$temp_output_file" > "$fasta_file"
        echo "完成: $fasta_file"
    else
        echo "警告: 处理 $fasta_file 后生成了空文件，跳过覆盖。"
    fi
    
    # 删除临时文件
    rm "$temp_output_file"
    
    echo "--------------------------------"
done

echo "所有 FASTA 文件处理完毕！"