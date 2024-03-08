from subprocess import run, PIPE
from datetime import datetime
import multiprocessing as mp
from send_email import send_email

FILE = "brand"
OUTPUT = "./" + FILE + ".out"
SOURCE = FILE + ".cpp"
COMPILE = ["g++", "-std=c++2a", "-I", "./eigen-3.4.0/", SOURCE, "-o", OUTPUT]

def function(i) :
    DIST = [.05, .75, 1., 1.25, 1.5]
    run(["mkdir", "-p", f"../data/{FILE}_D_{i}"], check=True)
    for r in DIST :
        with open("data_" + FILE + '_' + str(i) + '.txt', "w", encoding='utf-8') as data :
            data.write("100\n")
            data.write("2000000\n")
            data.write(f"../data/{FILE}_D_{i}/\n")
            data.write("40\n")
            data.write(f"{r}\n")

        if i == 0 :
            run([OUTPUT, str(i)], check=True)
        else :
            run([OUTPUT, str(i)], stdout=PIPE, stderr=PIPE, check=True)


if __name__ == "__main__" :
  start = datetime.now()
  run(COMPILE, check=True)

  pool = mp.Pool(20)
  pool.map(function, range(20))
  pool.join()

  end = datetime.now()
  text = "L'esecuzione_del_programma_" + FILE + "_D_e'_avvenuta_con_successo"
  send_email(subject="Fine esecuzione", body=text)
  print(text, "\nÈ iniziato alle ", start.strftime('%H:%M:%S'), "\ne si è concluso alle ", end.strftime('%H:%M:%S'), ".\n")
