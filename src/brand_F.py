from subprocess import run
from datetime import datetime
from send_email import send_email

FILE = "brand"
start = datetime.now()
OUTPUT = "./" + FILE + ".out"
SOURCE = FILE + ".cpp"
COMPILE = ["g++", "-std=c++2a", "-Wall", "-Wextra", "-Wpedantic", "-I", "./eigen-3.4.0/", SOURCE, "-o", OUTPUT]

run(COMPILE, check=True)
FORCES = [0, 20, 40, 60, 80, 100]
for i in range(0, 20) :
    run(["mkdir", f"./{FILE}_F_{i}"], check=True)
    for r in FORCES :
        with open("data_" + FILE + '.txt', "w", encoding='utf-8') as data :
            data.write("100\n")
            data.write("2000000\n")
            data.write(f"./{FILE}_F_{i}/\n")
            data.write(f"{r}\n")
            data.write("0.7\n")

        run([OUTPUT], check=True)

end = datetime.now()
text = "L'esecuzione_del_programma_" + FILE + "_F_e'_avvenuta_con_successo"
send_email("Fine esecuzione", text)
print(text, "\nÈ iniziato alle ", start.strftime('%H:%M:%S'), "\ne si è concluso alle ", end.strftime('%H:%M:%S'), ".\n")
