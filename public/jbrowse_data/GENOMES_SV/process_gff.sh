#!/bin/bash

# --- 配置 ---
# TSV 数据库文件，包含所有物种的 feature 信息
TSV_DATABASE="/home/data/CYM/database/43_species_features.tsv"
# 包含所有待处理 GFF 文件的根目录
GFF_DIR="/home/data/CYM/Webserver/xMetaVar/src/public/jbrowse_data/GENOMES_SV"

# --- 检查文件和目录是否存在 ---
if [ ! -f "$TSV_DATABASE" ]; then
    echo "错误：数据库文件 $TSV_DATABASE 未找到！"
    exit 1
fi

if [ ! -d "$GFF_DIR" ]; then
    echo "错误：GFF 目录 $GFF_DIR 未找到！"
    exit 1
fi

# --- 主逻辑 ---
echo "开始处理 GFF 文件..."

# 查找 GFF_DIR 下所有的 .gff 文件
find "$GFF_DIR" -type f -name "*.PATRIC.gff" | while IFS= read -r gff_file; do
    echo "正在处理: $gff_file"
    
    # 定义临时输出文件名
    temp_output_file="${gff_file}.tmp"
    
    # 核心 AWK 命令
    # -F'\t'       : 设置输入字段分隔符为 Tab
    # -v OFS='\t'    : 设置输出字段分隔符为 Tab
    # '...'          : awk 脚本主体
    # $TSV_DATABASE : 第一个输入文件 (查找表)
    # $gff_file     : 第二个输入文件 (待处理文件)
    awk -F'\t' -v OFS='\t' '
        # NR==FNR 这个条件只在读取第一个文件 (TSV_DATABASE) 时成立
        # NR 是总行号, FNR 是当前文件行号
        NR==FNR {
            # 跳过 TSV 文件的表头
            if (NR > 1) {
                # 使用 patric_id (第5列) 作为 key
                # 将 genome_id (第1列), global_start (第11列), global_end (第12列) 存入数组
                # db 数组存储了所有需要替换的信息
                patric_id = $5
                genome_id = $1
                global_start = $11
                global_end = $12
                db[patric_id, "gid"] = genome_id
                db[patric_id, "start"] = global_start
                db[patric_id, "end"] = global_end
            }
            # next 表示处理完当前行后，直接跳到下一行，不再执行后面的代码块
            next
        }

        # 这部分代码只在处理第二个文件 (gff_file) 时执行
        {
            # 如果行以 # 开头，是注释行，直接打印，不作处理
            if ($0 ~ /^#/) {
                print
                next
            }

            # 对于数据行，从第9列的属性中提取 ID
            # match函数会在$9中查找正则表达式，并将匹配结果存入 RSTART 和 RLENGTH
            # 我们用 substr 来提取匹配到的 ID=... 部分
            # 进一步用 sub 去掉 "ID=" 前缀
            if (match($9, /ID=[^;]+/)) {
                feature_id = substr($9, RSTART, RLENGTH)
                sub(/ID=/, "", feature_id)

                # 检查这个 ID 是否在我们的数据库 (db数组) 中
                if ((feature_id, "gid") in db) {
                    # 如果找到了，就执行替换操作
                    $1 = db[feature_id, "gid"]      # 替换第1列为 genome_id
                    $4 = db[feature_id, "start"]    # 替换第4列为 global_start
                    $5 = db[feature_id, "end"]      # 替换第5列为 global_end
                }
            }
            
            # 打印处理后 (或未处理) 的行
            print
        }
    ' "$TSV_DATABASE" "$gff_file" > "$temp_output_file"

    # 用处理后的临时文件覆盖原始文件
    mv "$temp_output_file" "$gff_file"
    echo "完成: $gff_file"
    echo "--------------------------------"
done

echo "所有 GFF 文件处理完毕！"