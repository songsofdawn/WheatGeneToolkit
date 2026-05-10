import sqlite3

# 注意：这是旧版单库测试脚本，当前主程序已使用 data/db 分库结构。
# 保留它仅用于检查历史 wheat_toolkit.db，不代表线上 Streamlit 的查询路径。
DB_FILE = "wheat_toolkit.db"

def main():
    conn = sqlite3.connect(DB_FILE)
    cur = conn.cursor()

    # 1. 查看所有表
    print("=== 数据库中的表 ===")
    cur.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;")
    tables = cur.fetchall()

    if not tables:
        print("数据库中没有表")
        return

    for t in tables:
        print(t[0])

    # 2. 依次查看每个表的结构和前5行
    for t in tables:
        table_name = t[0]
        print(f"\n\n=== 表名: {table_name} ===")

        # 查看表结构
        print("[表结构]")
        cur.execute(f"PRAGMA table_info({table_name})")
        cols = cur.fetchall()
        for col in cols:
            cid, name, col_type, notnull, dflt_value, pk = col
            print(f"  {name} ({col_type})")

        # 查看记录数
        cur.execute(f"SELECT COUNT(*) FROM {table_name}")
        n = cur.fetchone()[0]
        print(f"[记录数] {n}")

        # 查看前5行
        print("[前5行数据]")
        cur.execute(f"SELECT * FROM {table_name} LIMIT 5")
        rows = cur.fetchall()
        for row in rows:
            print(row)

    conn.close()

if __name__ == "__main__":
    main()
