from subprocess import run, PIPE
from datetime import datetime
import multiprocessing as mp
from send_email import send_email

FILE = "rod_1"
OUTPUT = "./" + FILE + ".out"
SOURCE = FILE + ".cpp"
COMPILE = ["g++", "-std=c++2a", "-Wall", "-Wextra", "-Wpedantic", "-I", "./eigen-3.4.0/", SOURCE, "-o", OUTPUT]

def function(i) :
    RADIUS = [2.1, 3.1, 4.1]
    run(["mkdir", f"./{FILE}_{i}"], check=True)
    for r in RADIUS :
        with open("data_" + FILE + '_' + str(i) + '.txt', "w", encoding='utf-8') as data :
            data.write("100\n")
            data.write("2000000\n")
            data.write(f"./{FILE}_{i}/\n")
            data.write(f"{r}\n")

        run([OUTPUT, str(i)], stdout=PIPE, check=True)

if __name__ == "__main__" :
    start = datetime.now()
    run(COMPILE, check=True)
    
    pool = mp.Pool(20)
    pool.map(function, range(20))
    pool.close()

    end = datetime.now()
    text = "L'esecuzione_del_programma_" + FILE + "_e'_avvenuta_con_successo"
    send_email("Fine esecuzione", text)
    print(text, "\nÈ iniziato alle ", start.strftime('%H:%M:%S'), "\ne si è concluso alle ", end.strftime('%H:%M:%S'), ".\n")
