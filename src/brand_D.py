from subprocess import run
from datetime import datetime
from send_email import send_email

FILE = "brand"
start = datetime.now()
OUTPUT = "./" + FILE + ".out"
SOURCE = FILE + ".cpp"
COMPILE = ["g++", "-std=c++2a", "-Wall", "-Wextra", "-Wpedantic", "-I", "./eigen-3.4.0/", SOURCE, "-o", OUTPUT]

# run(COMPILE, check=True)
# DIST = [.05, .75, 1., 1.25, 1.5]
# for i in range(0, 20) :
#     run(["mkdir", f"./{FILE}_D_{i}"], check=True)
#     for r in DIST :
#         with open("data_" + FILE + '.txt', "w", encoding='utf-8') as data :
#             data.write("100\n")
#             data.write("2000000\n")
#             data.write(f"./{FILE}_D_{i}/\n")
#             data.write("40\n")
#             data.write(f"{r}\n")

#         run([OUTPUT], check=True)

end = datetime.now()
text = "L'esecuzione_del_programma_" + FILE + "_D_e'_avvenuta_con_successo"
send_email(subject="Fine esecuzione", body=text)
print(text, "\nÈ iniziato alle ", start.strftime('%H:%M:%S'), "\ne si è concluso alle ", end.strftime('%H:%M:%S'), ".\n")
