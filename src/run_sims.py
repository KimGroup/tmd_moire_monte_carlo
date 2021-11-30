from datetime import date
from string import ascii_lowercase
import os

# function to run several simulations in the background across multiple cores
def main():
  mode = "new"
  prefix = date.today().strftime("%d%m%y") # output file prefix
  f = open("runs.log",'r')
  lines = f.readlines()
  f.close()
  exists = []
  # add a letter to make the prefix unique if multiple runs on same date
  for line in lines:
    if prefix in line: exists += [line]
  if len(exists) > 0:
    if mode == "append":
      raise NotImplementedError()
    else:
      for letter in ascii_lowercase:
        if len([i for i in filter(lambda x: prefix+letter in x,exists)]) == 0:
          prefix += letter
          break
  else: prefix += 'a'
  # log the prefix, then go into runs.log to add notes/details
  f = open("runs.log",'a')
  f.writelines(["\nRUN: "+prefix+"\n"])
  f.close()
  
  # lattice size and list of filling fraction numerators to run

  #l,nums = 4,[6, 7, 8]
  #l,nums = 12,[48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72]
  l, nums = 9,[36]
  #l,nums = 10,[34,36,38,40,42,44,46,48,50]
  #l,nums = 20,[134, 136, 144, 152, 160, 168, 176, 184, 192, 200]
  #l, nums = 30,[300,306,324,360,378,396,414,432,450]
  #l,nums = 6,[12,13,14,15,16,17,18]
  
  # submit the simulations
  command = "nohup python particle_mc.py {l} {num} {pref} &"
  for n in nums: os.system(command.format(l=l,num=n,pref=prefix))

if __name__ == '__main__': main()
