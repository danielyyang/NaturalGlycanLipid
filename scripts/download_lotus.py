"""
LOTUS 数据库 Dump 下载器 — 支持代理与断点续传
LOTUS Database Dump Downloader — Proxy support with retry mechanism

下载 Zenodo 上的 LOTUS frozen metadata CSV (230106_frozen_metadata.csv.gz, ~108MB)
"""
import os
import sys
import time
import requests
from tqdm import tqdm
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# =====================================================================
# ⚙️ 配置区 — 请根据你的代理环境修改以下参数
# ⚙️ Configuration — Modify these parameters for your proxy setup
# =====================================================================

# 代理地址 (Proxy address) — 设为 None 可禁用代理
PROXY_HOST = "127.0.0.1"
PROXY_PORT = 7890

# 目标 URL (Target URL)
LOTUS_URL = "https://zenodo.org/records/7534071/files/230106_frozen_metadata.csv.gz"

# 最大重试次数 (Max retries)
MAX_RETRIES = 3

# 连接超时 / 读取超时 (秒)
CONNECT_TIMEOUT = 30
READ_TIMEOUT = 120

# 下载分块大小 (Download chunk size in bytes)
CHUNK_SIZE = 8192


def buildSession(proxyHost: str = PROXY_HOST, proxyPort: int = PROXY_PORT) -> requests.Session:
    """
    构建带代理和自动重试的 requests Session。
    Build a requests Session with proxy and automatic retry support.
    """
    session = requests.Session()

    # 配置重试策略 (Configure retry strategy)
    retryStrategy = Retry(
        total=MAX_RETRIES,
        backoff_factor=1,            # 重试间隔: 1s, 2s, 4s (Retry backoff: 1s, 2s, 4s)
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "HEAD"],
    )
    adapter = HTTPAdapter(max_retries=retryStrategy)
    session.mount("https://", adapter)
    session.mount("http://", adapter)

    # 配置代理 (Configure proxy)
    if proxyHost and proxyPort:
        proxyUrl = f"http://{proxyHost}:{proxyPort}"
        session.proxies = {
            "http": proxyUrl,
            "https": proxyUrl,
        }
        print(f"  Proxy: {proxyUrl}")
    else:
        print("  Proxy: disabled (direct connection)")

    return session


def downloadWithProgress(
    url: str,
    outputPath: str,
    session: requests.Session,
) -> bool:
    """
    带进度条和重试机制的文件下载器。
    File downloader with tqdm progress bar and retry mechanism.

    Args:
        url: 下载 URL
        outputPath: 保存路径
        session: 配置好的 requests Session

    Returns:
        bool: 下载是否成功
    """
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            print(f"\n  Attempt {attempt}/{MAX_RETRIES}: Connecting to {url}...")
            response = session.get(
                url,
                stream=True,
                timeout=(CONNECT_TIMEOUT, READ_TIMEOUT),
            )
            response.raise_for_status()

            # 获取文件总大小 (Get total file size)
            totalSize = int(response.headers.get("content-length", 0))
            totalSizeMb = totalSize / (1024 * 1024)
            print(f"  File size: {totalSizeMb:.1f} MB")

            # 下载并写入文件 (Download and write to file)
            downloadedSize = 0
            with open(outputPath, "wb") as f:
                with tqdm(
                    total=totalSize,
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
                    desc=f"  Downloading",
                    ncols=80,
                ) as progressBar:
                    for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
                        if chunk:
                            f.write(chunk)
                            downloadedSize += len(chunk)
                            progressBar.update(len(chunk))

            # 验证完整性 (Verify completeness)
            if totalSize > 0 and downloadedSize != totalSize:
                print(f"  WARNING: Downloaded {downloadedSize} bytes, expected {totalSize} bytes.")
                print(f"  File may be incomplete. Retrying...")
                continue

            print(f"  Download complete: {outputPath}")
            return True

        except requests.exceptions.ConnectionError as e:
            print(f"  Connection error on attempt {attempt}: {e}")
            if attempt < MAX_RETRIES:
                waitTime = 2 ** attempt
                print(f"  Retrying in {waitTime}s...")
                time.sleep(waitTime)
        except requests.exceptions.Timeout as e:
            print(f"  Timeout on attempt {attempt}: {e}")
            if attempt < MAX_RETRIES:
                waitTime = 2 ** attempt
                print(f"  Retrying in {waitTime}s...")
                time.sleep(waitTime)
        except requests.exceptions.HTTPError as e:
            print(f"  HTTP error on attempt {attempt}: {e}")
            if attempt < MAX_RETRIES:
                waitTime = 2 ** attempt
                print(f"  Retrying in {waitTime}s...")
                time.sleep(waitTime)
        except KeyboardInterrupt:
            print("\n  Download cancelled by user.")
            # 清理不完整的文件 (Clean up incomplete file)
            if os.path.exists(outputPath):
                os.remove(outputPath)
                print(f"  Removed incomplete file: {outputPath}")
            return False

    print(f"  FAILED: All {MAX_RETRIES} attempts exhausted.")
    return False


def main():
    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    dataDir = os.path.join(baseDir, "data")
    os.makedirs(dataDir, exist_ok=True)

    outputPath = os.path.join(dataDir, "230106_frozen_metadata.csv.gz")

    # 检查是否已下载 (Check if already downloaded)
    if os.path.exists(outputPath):
        existingSize = os.path.getsize(outputPath) / (1024 * 1024)
        print(f"File already exists: {outputPath} ({existingSize:.1f} MB)")
        userInput = input("Re-download? [y/N]: ").strip().lower()
        if userInput != "y":
            print("Skipped. Use existing file.")
            return

    print("=" * 70)
    print("LOTUS Frozen Metadata Downloader")
    print("=" * 70)
    print(f"  URL: {LOTUS_URL}")
    print(f"  Output: {outputPath}")
    print(f"  Max retries: {MAX_RETRIES}")

    session = buildSession()
    success = downloadWithProgress(LOTUS_URL, outputPath, session)

    if success:
        finalSize = os.path.getsize(outputPath) / (1024 * 1024)
        print(f"\n{'=' * 70}")
        print(f"SUCCESS! Downloaded {finalSize:.1f} MB to {outputPath}")
        print(f"{'=' * 70}")
        print(f"\nNext step — run taxonomy filling:")
        print(f"  python -c \"from lib.taxonomy_lotus_matcher import runLocalTaxonomyFilling; "
              f"runLocalTaxonomyFilling('reports/Coconut_Sugar_Check.csv', "
              f"'data/230106_frozen_metadata.csv.gz', "
              f"'reports/Coconut_Sugar_WithTax.csv')\"")
    else:
        print(f"\nFAILED to download. Please check your proxy settings at the top of this script.")


if __name__ == "__main__":
    main()
