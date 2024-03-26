from subprocess import run, PIPE
from datetime import datetime
import multiprocessing as mp
from send_email import send_email

FILE = "histone"
OUTPUT = "./" + FILE + ".out"
SOURCE = FILE + ".cpp"
COMPILE = ["g++", "-std=c++2a", "-I", "./eigen-3.4.0/", SOURCE, "-o", OUTPUT]

def function(i) :
    run(["mkdir", "-p", f"../data/{FILE}_{i}"], check=True)
    with open("data_" + FILE + '_' + str(i) + '.txt', "w", encoding='utf-8') as data :
        data.write("200\n")
        data.write("600000\n")
        data.write(f"../data/{FILE}_{i}/\n")
        data.write("0.\t-4.5\t10.2\n")
    
    if i == 0 :
        run([OUTPUT, str(i)], check=True)
    else :
        run([OUTPUT, str(i)], stdout=PIPE, stderr=PIPE, check=True)

if __name__ == "__main__" :
    send_email("Inizio esecuzione", f"L'esecuzione_del_programma_{FILE}_e'_iniziata")
    start = datetime.now()
    run(COMPILE, check=True)

    pool = mp.Pool(10)
    pool.map(function, range(10))
    pool.close()
    pool.join()

    end = datetime.now()
    text = "L'esecuzione_del_programma_" + FILE + "_e'_avvenuta_con_successo"
    send_email("Fine esecuzione", text)
    print(text, "\nÈ iniziato alle ", start.strftime('%H:%M:%S'),
          "\ne si è concluso alle ", end.strftime('%H:%M:%S'), ".\n")
