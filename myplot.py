import pandas as pd
import matplotlib.pyplot as plt
import sys

# Komut satırından dosya adını alın
if len(sys.argv) != 2:
    print("Kullanım: python plot.py <dosya_adı>")
    sys.exit(1)

dosya_adı = sys.argv[1]

# Veriyi pandas DataFrame'e oku (ayrıştırma için boşlukları kullan)
veri = pd.read_table(dosya_adı, sep="\s+", header=None, 
       names=["iter","reward" "reward50"],index_col="iter",usecols=[0,1])

plt.plot(veri)

# Grafiği özelleştirme
plt.title("Veri Grafiği")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()

# Grafiği gösterme
plt.show()
