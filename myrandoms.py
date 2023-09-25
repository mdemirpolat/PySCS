import random

# Rastgele sayıları üretip dosyaya yazma
# DOS'da kullanabilmek için Line Ending düzeltilmeli (in Sublime : View/LineEnding)



n = 5000  # Üreteceğiniz sayı adedi
filename = "RNDHAZIR.TXT"
randomfinish = False
orneksayisi=10

def generate_and_write_random_numbers(filename, n):
    with open(filename, "w") as file:
        for _ in range(n):
            random_num = random.uniform(0, 1)
            file.write(str(random_num) + "\n")
    print(f"(0-1) aralığında {n} tane random sayı üretildi ve {filename} dosyasına başarıyla yazıldı.")

# Dosyadan sırasıyla okuma ve rastgele sayı üreteci olarak kullanma
def frandom(filename):  
    while True and not randomfinish:
        with open(filename, "r") as file:
            lines = file.readlines()
            for line in lines:
                yield float(line.strip())

# Kullanım
generate_and_write_random_numbers(filename, n)
random_gen = frandom(filename)
k=0
print(f"{orneksayisi} kadar ornek asagidadir :")
for _ in range(orneksayisi):
    k+=1
    random_num = next(random_gen)
    print(k,random_num)

