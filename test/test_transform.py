import numpy as np
from lookup_tables import make_transformation_table, make_nn_table

def main():
  l = 20
  t_table = make_transformation_table(l)
  nn_table = make_nn_table((l,l),max_neighbor=1)
  for i in range(t_table.shape[0]):
    for j in range(l*l):
      if set(t_table[i,:][nn_table[j]]) != set(nn_table[t_table[i,j]]):
        print("Transformation failed for transformation {t} at site {s}".format(t=i,s=j))
  return

if __name__ == '__main__': main()
