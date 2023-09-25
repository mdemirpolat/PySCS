import pandas as pd
import matplotlib.pyplot as plt
import sys

# Komut satırından dosya adını alın
if len(sys.argv) != 2:
    print("Kullanım: python plot.py <dosya_adı>")
    sys.exit(1)

dosya_adı = sys.argv[1]

# Veriyi pandas DataFrame'e oku (ayrıştırma için boşlukları kullan)
veri = pd.read_csv(dosya_adı, delimiter="\s+")

# İlk sütunu x değeri olarak al
x = veri.iloc[:, 0]

# Diğer sütunları y değerleri olarak al
y_sutunlar = veri.columns[1:]

# Yeni bir sütun adları listesi oluştur
y_sutun_adlari = [f"y{i+1}" for i in range(len(y_sutunlar))]

# Veriyi kopyala ve yeni sütun adlarını kullanarak sütunları yeniden adlandır
veri_yeniden_adlandirilmis = veri.copy()
veri_yeniden_adlandirilmis.columns = ["x"] + y_sutun_adlari

# Her y sütununu çizdir
for i, y_sutun in enumerate(y_sutun_adlari):
    if i == 0:  # İlk sütun (y1) kalın çizgi
        plt.plot(x, veri_yeniden_adlandirilmis[y_sutun], label=y_sutun, linewidth=2)
    else:  # Diğer sütunlar (y2, y3, ...) noktalı ve daha ince çizgi
        plt.plot(x, veri_yeniden_adlandirilmis[y_sutun], label=y_sutun, linestyle='--', alpha=0.5)

# Grafiği özelleştirme
plt.title("Veri Grafiği")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()

# Grafiği gösterme
plt.show()
