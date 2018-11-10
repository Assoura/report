import os
def push(wait=300,path = os.path.join('/','home','pi','Bureau','Mirmoc_is_back')):
    os.chdir(path)
    run=['git add temp.png','git commit -m "maj"','git push -f']
    for c in run: os.system(c)

push()