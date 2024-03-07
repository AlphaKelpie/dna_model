from subprocess import run
from datetime import datetime
from send_email import send_email

FILE = "histone"
start = datetime.now()
OUTPUT = "./" + FILE + ".out"
SOURCE = FILE + ".cpp"
COMPILE = ["g++", "-std=c++2a", "-Wall", "-Wextra", "-Wpedantic", "-I", "./eigen-3.4.0/", SOURCE, "-o", OUTPUT]

run(COMPILE, check=True)

for i in range(0, 10) :
    with open("data_" + FILE + '.txt', "w", encoding='utf-8') as data :
        data.write("200\n")
        data.write("600000\n")
        data.write(f"./{FILE}_{i}/\n")
        data.write("0.\t-4.5\t10.2\n")
    
    run(["mkdir", f"./{FILE}_{i}"], check=True)
    run([OUTPUT], check=True)

end = datetime.now()
text = "L'esecuzione_del_programma_" + FILE + "_e'_avvenuta_con_successo"
send_email("Fine esecuzione", text)
print(text, "\nÈ iniziato alle ", start.strftime('%H:%M:%S'), "\ne si è concluso alle ", end.strftime('%H:%M:%S'), ".\n")
