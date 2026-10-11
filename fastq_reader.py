# -*- coding: utf-8 -*-
"""独立 FASTQ 读取器：逐条返回完整记录，只依赖 Python 标准库。"""

# 【名称与作用】
# - 名称：read_fastq；本文件可独立上传 GitHub、复制使用，不依赖 run.py 或教学文件。
# - 作用：准确读取记录的标题、序列、+ 行和原始质量字符串；已确认编码时可转换逐碱基 Phred 分值。
# - 流程位置：已有 FASTQ → 逐条解析及结构检查 → 由调用者使用序列和质量开展后续分析。
#
# 【状态与验证】
# - 代码状态：已验证（2026-10-10，qa/test_standalone_fastq.py 的 14 项模拟回归全部通过）。
# - 运行情况与正确性核对：本次实际运行，用人工预期序列、标题、质量字符及分值逐项比较；
#   将此文件单独复制到不含 run.py 的目录，用 Python 隔离模式执行成功。
# - 验证范围：普通/gzip、可变读长、折行、Phred+33/+64、原文模式、坏记录、gzip 完整性、
#   部分读取及 1,000/20,000 条模拟记录的内存检查；大型真实数据、跨平台和独立解析器对照未验证。
#
# 【适用条件与限制】
# - 适用条件：普通/gzip 的 ASCII DNA FASTQ，可变读长；可配置接受四行式或折行式记录。
# - 不支持或未检查：不解码 Solexa 分值、颜色空间、RNA 的 U、其他压缩格式或 BAM/CRAM；
#   不检查 R1/R2 配对、read 名称唯一性、接头或污染，不做过滤、比对及变异调用。
# - 不能推出：记录可读不代表完整质控通过；仅看前 N 条不代表完整文件，Phred 字符合法不证明编码选对。
#
# 【环境与依赖】
# - 语言和包：仅用 Python 标准库 pathlib、gzip、zlib、json、contextlib。
# - 版本要求：Python 3.8 或更新版本作为使用范围；最低版本尚未实测。
# - 实际验证版本、系统和架构：Python 3.13.7、macOS 27.0.1、arm64；来自本次测试进程记录。
# - 参考和上游变量：不需要参考基因组、注释库或其他项目文件；输入来自 path 或下方 INPUT_FASTQ。
#
# 【输入与参数】
# - 输入：字符串/Path 路径；记录必须包含非空 @ 标题、非空序列、+ 行和等长质量字符串。
# - 需确认：文件路径、允许的碱基与折行规则；phred_offset=None 只保留原始质量，不猜测编码；
#   确认是 Phred+33/+64 时填 33/64，才计算数值质量。None 不等于已识别为某种质量编码。
# - 路径：推荐绝对路径，支持 ~；相对路径按首次迭代时的工作目录解析，不依赖 __file__。
#
# 【输出与影响】
# - 输出对象：生成器逐条产生字典，字段是 title、sequence、plus、quality、phred_scores。
# - 字段：title 为 @ 后的完整标题；sequence 保留碱基原有大小写；plus 保留完整 + 行；
#   quality 保留原始质量字符；phred_scores 为等长整数列表，未指定编码时为 None。
# - 折行：合并同一记录的序列/质量行，只去掉行结束符；不保留原文件折行排版，不是逐字节文件副本。
# - 副作用：只读输入，不创建或覆盖数据文件；直接运行会打印少量记录预览和读取范围；
#   导入无数据读写副作用。生成器只能消费一次，提前退出请使用 closing() 或显式 close()。
#
# 【运行与检查】
# - 入口：填写 INPUT_FASTQ 后完整运行本文件；其他脚本可 from fastq_reader import read_fastq。
# - 最小调用：for record in read_fastq("/待替换/sample.fastq.gz", phred_offset=None): ...
#   路径需替换，循环体需填写自己的处理逻辑；不要用 list(...) 把大型文件整体放进内存。
# - 优先检查：文件路径、读取范围、首条记录的序列和质量长度、编码是否有来源依据。
# - 预期：例如序列 AC、质量 I!，在确认 Phred+33 时对应 [40, 0]；未指定编码时仍保留 I!，数值为 None。
# - 常见错误：文件不存在查路径，格式错误按记录/行号定位；不要随意放宽编码和碱基规则绕过坏数据。
#   max_reads=None 也须完整耗尽生成器且没有异常，才算读到末尾；提前 break 不能验证文件尾部。
#
# 【依据与待确认事项】
# - 格式依据：https://biopython.org/docs/latest/api/Bio.SeqIO.QualityIO.html
# - gzip 读取：https://docs.python.org/3.13/library/gzip.html （实现未使用 Biopython）。
# - 待确认：实际样本格式、编码、配对状态及大文件性能；本工具只承诺在已声明且已测试的格式范围内读取。

from pathlib import Path
from contextlib import closing
import gzip
import json
import zlib


# ==================== 使用时先检查这个配置区 ====================
INPUT_FASTQ = None       # 【必填】例如 "/absolute/path/sample_R1.fastq.gz"，替换为实际路径。
PHRED_OFFSET = None      # 【需确认】None=仅保留质量原文；确认 Phred+33/+64 后填 33/64，不支持 Solexa 数值转换。
MAX_READS = 10000         # 【需确认范围】预览前一万条；检查完整文件设 None，前 N 条不是随机抽样。
ALLOW_WRAPPED = False    # 【需确认格式】False=四行式；确实有序列/质量折行时设 True。
ALLOWED_BASES = "ACGTN"   # 【需确认字符】若来源包含 DNA IUPAC 模糊碱基，可填 "ACGTRYSWKMBDHVN"。
COMPRESSION = "auto"     # 【通常保留】按文件头识别 gzip；也可明确指定 "plain" 或 "gzip"。
ALLOW_EMPTY = False      # 【通常保留】空文件报错；只有预期允许空输入时设 True。
PREVIEW_RECORDS = 3      # 【显示设置】完整打印前几条记录，0=不打印记录；不改变 MAX_READS 的读取范围。
# 本文件只负责读取，不再提供 CSV 写出/覆盖参数。返回的完整记录可交给你后续的处理代码。


class FastqFormatError(ValueError):
    """输入不符合声明的 FASTQ 格式时，指出文件及出错位置。"""


def _compression_type(path, compression):
    """按文件头确认压缩方式，不根据扩展名猜测文件内容。"""
    if compression not in ("auto", "plain", "gzip"):
        raise ValueError("compression 只能是 'auto'、'plain' 或 'gzip'")
    with path.open("rb") as handle:
        detected = "gzip" if handle.read(2) == b"\x1f\x8b" else "plain"
    if compression != "auto" and compression != detected:
        raise ValueError(f"{path}：文件头与 compression={compression!r} 不符")
    return detected


def read_fastq(
    path, *, phred_offset=None, compression="auto", allow_wrapped=False,
    allowed_bases="ACGTN", max_reads=None, allow_empty=False,
):
    """逐条返回完整 FASTQ 记录，不累积整个文件，不静默跳过损坏记录。

    每条返回：
        {"title": "read1", "sequence": "AC", "plus": "+",
         "quality": "I!", "phred_scores": [40, 0]}
    上例的分值只适用于 phred_offset=33；None 时 phred_scores 为 None。

    path：必填输入路径。其余参数需用名称填写，含义见顶部配置区。
    本函数是生成器，首次迭代时才检查参数和打开输入；只构造生成器并不会读取文件。
    内存主要取决于单条 read 大小；数值质量启用时也只为当前 read 建立分值列表。
    """
    path = Path(path).expanduser().resolve()
    if phred_offset is not None and (type(phred_offset) is not int or phred_offset not in (33, 64)):
        raise ValueError("phred_offset 只能是 None、33 或 64；None 表示不解码质量分值")
    if max_reads is not None and (type(max_reads) is not int or max_reads <= 0):
        raise ValueError("max_reads 必须是正整数或 None")
    if not isinstance(allow_wrapped, bool) or not isinstance(allow_empty, bool):
        raise ValueError("allow_wrapped 和 allow_empty 必须为 True 或 False")
    if (not isinstance(allowed_bases, str) or not allowed_bases
            or any(c.upper() not in "ACGTRYSWKMBDHVN" for c in allowed_bases)):
        raise ValueError("allowed_bases 需为非空 DNA/IUPAC 字母字符串，例如 'ACGTN'")
    bases = set(allowed_bases.upper())
    opener = gzip.open if _compression_type(path, compression) == "gzip" else open
    record_number, line_number = 0, 0

    def fail(message):
        return FastqFormatError(
            f"{path}：第 {record_number} 条 read，附近第 {line_number} 行：{message}"
        )

    try:
        with opener(path, "rt", encoding="ascii", newline=None) as handle:
            def take_line():
                nonlocal line_number
                text = handle.readline()  # 只读一行；空行与文件结束必须区分。
                if text == "":
                    return None
                line_number += 1
                return text.rstrip("\r\n")  # 只去换行，不使用 strip() 隐藏非法空格。

            while max_reads is None or record_number < max_reads:
                header = take_line()
                if header is None:
                    break
                record_number += 1
                if not header.startswith("@") or not header[1:].strip():
                    raise fail("名称行必须以 @ 开头，且标题不能为空")
                title = header[1:]  # 不拆分空格后的说明，保留完整标题。

                if allow_wrapped:
                    sequence_parts = []
                    while True:
                        line = take_line()
                        if line is None:
                            raise fail("文件被截断，缺少 + 分隔行")
                        if line.startswith("+"):
                            plus = line
                            break
                        if not line or any(c.upper() not in bases for c in line):
                            raise fail("序列含空行或不在 allowed_bases 中的字符")
                        sequence_parts.append(line)
                    sequence = "".join(sequence_parts)
                else:
                    sequence, plus = take_line(), take_line()
                    if sequence is None or plus is None:
                        raise fail("文件被截断，缺少序列或 + 分隔行")
                    if not plus.startswith("+"):
                        raise fail("预期为 + 分隔行；请检查文件，或确认是否确需 allow_wrapped=True")

                if not sequence:
                    raise fail("不接受空序列")
                invalid = sorted(set(sequence.upper()) - bases)
                if invalid:
                    raise fail(f"序列含不允许的字符 {invalid!r}；请核对来源和 allowed_bases")
                if plus[1:] and plus[1:] != title:
                    raise fail("+ 后重复的标题与 @ 行标题不一致")

                quality_parts, quality_length = [], 0
                # 未确认编码时仍要求质量为可打印的非空白 ASCII 字符，不自行解释其数值。
                minimum_code = 33 if phred_offset is None else phred_offset
                while quality_length < len(sequence):
                    line = take_line()
                    if line is None:
                        raise fail("文件被截断，质量字符串不足")
                    if not line:
                        raise fail("质量字符串为空")
                    for char in line:
                        if not minimum_code <= ord(char) <= 126:
                            raise fail(f"质量字符 {char!r} 不符合声明的编码/ASCII 范围")
                    quality_parts.append(line)
                    quality_length += len(line)
                    if not allow_wrapped:
                        break
                if quality_length != len(sequence):
                    raise fail(f"序列长度 {len(sequence)} 与质量长度 {quality_length} 不一致")
                quality = "".join(quality_parts)
                scores = None if phred_offset is None else [ord(c) - phred_offset for c in quality]

                # 每次交出一条完整记录；原始字符不因质量解码而被替换或修改。
                yield {
                    "title": title, "sequence": sequence, "plus": plus,
                    "quality": quality, "phred_scores": scores,
                }
    except (OSError, EOFError, UnicodeError, zlib.error) as exc:
        raise FastqFormatError(
            f"{path}：第 {record_number} 条 read 附近读取失败（缓冲可能使行号不精确）：{exc}"
        ) from exc
    if record_number == 0 and not allow_empty:
        raise FastqFormatError(f"{path}：没有 read；只有明确允许空输入时才设 allow_empty=True")


def main():
    """直接运行时展示读取效果；在其他脚本导入本文件时不会执行。"""
    if INPUT_FASTQ is None:
        raise SystemExit("请先在配置区填写 INPUT_FASTQ；PHRED_OFFSET 未确认可保留 None，仅读取质量原文。")
    if type(PREVIEW_RECORDS) is not int or PREVIEW_RECORDS < 0:
        raise SystemExit("PREVIEW_RECORDS 必须是非负整数")
    count, total_bases = 0, 0
    records = read_fastq(
        INPUT_FASTQ, phred_offset=PHRED_OFFSET, compression=COMPRESSION,
        allow_wrapped=ALLOW_WRAPPED, allowed_bases=ALLOWED_BASES,
        max_reads=MAX_READS, allow_empty=ALLOW_EMPTY,
    )
    with closing(records):  # 即使提前退出或后续处理报错，也关闭已打开的输入文件。
        for record in records:
            count += 1
            total_bases += len(record["sequence"])
            if count <= PREVIEW_RECORDS:
                print(json.dumps({"type": "record_preview", "record": record}, ensure_ascii=False))
            # 【接入自己的处理】可在这里使用 record["sequence"]、record["quality"] 等字段。
            # 此循环不把全部记录存入列表；打印的记录预览不代表最后一定读取成功。
    print(json.dumps({
        "type": "completed", "input": str(Path(INPUT_FASTQ).expanduser().resolve()),
        "read_count": count, "total_bases": total_bases,
        "scan_mode": "FULL_FILE" if MAX_READS is None else "FIRST_N_READS",
        "whole_file_checked": MAX_READS is None,
        "phred_offset": PHRED_OFFSET, "quality_decoded": PHRED_OFFSET is not None,
        "max_reads": MAX_READS, "allow_wrapped": ALLOW_WRAPPED,
        "allowed_bases": ALLOWED_BASES, "compression_setting": COMPRESSION,
        "paired_end_checked": False,
    }, ensure_ascii=False, indent=2))
    # 只有输出 completed 且进程正常结束，才表示本次指定范围读取成功。
    # 设置 MAX_READS 后，即使文件恰好较短，也保守标记为 FIRST_N_READS。


if __name__ == "__main__":
    main()
