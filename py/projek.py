from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

# ==========================
# Membaca file FASTA
# ==========================

file_fasta = "hepartitis_b.fasta"

# Menyimpan data dalam List
data = []

for record in SeqIO.parse(file_fasta, "fasta"):

    seq = str(record.seq).upper()

    panjang = len(seq)

    # Dictionary frekuensi nukleotida
    frekuensi = {
        "A": seq.count("A"),
        "T": seq.count("T"),
        "G": seq.count("G"),
        "C": seq.count("C")
    }

    # Hitung GC Content
    gc_content = (
        (frekuensi["G"] + frekuensi["C"])
        / panjang
    ) * 100

    data.append({
        "ID": record.id,
        "Length": panjang,
        "A": frekuensi["A"],
        "T": frekuensi["T"],
        "G": frekuensi["G"],
        "C": frekuensi["C"],
        "GC_Content": round(gc_content, 2)
    })

# ==========================
# DataFrame
# ==========================

df = pd.DataFrame(data)

# Urutkan berdasarkan GC Content
df_sorted = df.sort_values(
    by="GC_Content",
    ascending=False
)

# ==========================
# Menampilkan jumlah sekuens
# ==========================

print(f"\nJumlah sekuens: {len(df_sorted)}")

# ==========================
# Top 3 Sekuens Terbaik
# ==========================

print("\n=== TOP 3 GC CONTENT ===")

top3 = df_sorted.head(3)

print(
    top3[["ID", "GC_Content"]]
)

# ==========================
# Simpan ke CSV
# ==========================

df_sorted.to_csv(
    "hasil.csv",
    index=False
)

# ==========================
# Visualisasi Grafik
# ==========================

plt.figure(figsize=(8,5))

plt.bar(
    df_sorted["ID"],
    df_sorted["GC_Content"]
)

plt.title("Perbandingan GC Content")
plt.xlabel("ID Sekuens")
plt.ylabel("GC Content (%)")

plt.xticks(rotation=30)

plt.tight_layout()

plt.savefig("hasil.png")

plt.show()

# ==========================
# Data Lengkap
# ==========================

print("\n=== DATA LENGKAP ===")
print(df_sorted)

print("\nFile CSV berhasil dibuat: hasil.csv")
print("Grafik berhasil dibuat: hasil.png")