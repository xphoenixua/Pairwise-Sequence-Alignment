import numpy as np
import pandas as pd

# для запам’ятовування шляху потрапляння в точку
def trace_M(i, j, idx):
    if idx==0:
        trace_point = ['M', (i-1, j-1)]
    elif idx==1:
        trace_point = ['Bx', (i-1, j-1)]
    elif idx==2:
        trace_point = ['By', (i-1, j-1)]
    return trace_point

def trace_Bx(i, j, idx):
    if idx==0:
        trace_point = ['M', (i-1, j)]
    elif idx==1:
        trace_point = ['Bx', (i-1, j)]
    return trace_point

def trace_By(i, j, idx):
    if idx==0:
        trace_point = ['M', (i, j-1)]
    elif idx==1:
        trace_point = ['By', (i, j-1)]
    return trace_point

# умови задачі
d = 4
e = 1
seq_1 = 'LGLI'
seq_2 = 'IQGL'
seq_1 = [0] + [*seq_1]
seq_2 = [0] + [*seq_2]

# ініціалізація матриці ваг та матриці зворотнього проходу
n, m = len(seq_1), len(seq_2)
M, Bx, By = np.zeros((n,m)), np.zeros((n,m)), np.zeros((n,m))
M_T = np.zeros((n,m), dtype=list)
Bx_T = np.zeros((n,m), dtype=list)
By_T = np.zeros((n,m), dtype=list)

# задання початкових умов матриць прямого та зворотнього проходу
M[0,0], Bx[0,0], By[0,0] = 0, -np.inf, -np.inf
Bx_T[1,0], By_T[0,1] = ['M', (0,0)], ['M', (0,0)]
for i in range(1, n):
    M[i,0], By[i,0] = -np.inf, -np.inf
    Bx[i,0] = -d - (i - 1)*e
    if i>1:
        Bx_T[i,0] = ['Bx', (i-1, 0)]
for j in range(1, m):
    M[0,j], Bx[0,j] = -np.inf, -np.inf
    By[0,j] = -d - (j - 1)*e
    if j>1:
        By_T[0,j] = ['By', (0, j-1)]

# заповнення комірок матриць
for j in range(1, m):
    for i in range(1, n):
        if seq_1[i]==seq_2[j]:
            s = 1
        else:
            s = -1
        max_M_arr = [M[i-1, j-1] + s, Bx[i-1, j-1] + s, By[i-1, j-1] + s]
        max_Bx_arr = [M[i-1, j] - d, Bx[i-1, j] - e]
        max_By_arr = [M[i, j-1] - d, By[i, j-1] - e]
        M[i,j] = np.max(max_M_arr)
        Bx[i,j] = np.max(max_Bx_arr)
        By[i,j] = np.max(max_By_arr)

        # залежно від варіанту запам’ятовую звідки було здійснено переміщення
        idx_M = np.argmax(max_M_arr)
        idx_Bx = np.argmax(max_Bx_arr)
        idx_By = np.argmax(max_By_arr)
        M_T[i,j] = trace_M(i, j, idx_M)
        Bx_T[i,j] = trace_Bx(i, j, idx_Bx)
        By_T[i,j] = trace_By(i, j, idx_By)

# відстеження зворотнього шляху від комірки F(n,m) до комірки F(0,0)
main_length = max([n,m])
trace = np.zeros(main_length, dtype=list)

# задаємо початкові умови для шляху зворотнього проходу
max_arr = [M[n-1,m-1], Bx[n-1,m-1], By[n-1,m-1]]
idx = np.argmax(max_arr)
if idx==0:
    trace[0] = ['M', (n-1,m-1)]
elif idx==1:
    trace[0] = ['Bx', (n-1,m-1)]
elif idx==2:
    trace[0] = ['By', (n-1,m-1)]

# відновлення зворотнього шляху
i, j = n-1, m-1
xy = 1
while i!=0 and j!=0:
    prev = trace[xy-1]
    if prev[0]=='M':
        trace[xy] = M_T[prev[1][0], prev[1][1]]
    elif prev[0]=='Bx':
        trace[xy] = Bx_T[prev[1][0], prev[1][1]]
    elif prev[0]=='By':
        trace[xy] = By_T[prev[1][0], prev[1][1]]
    i, j = trace[xy][1][0], trace[xy][1][1]
    xy += 1
trace = trace[::-1]

# вирівнювання двох послідовностей
seq_1_alignment = np.zeros(main_length, dtype=object)
seq_2_alignment = np.zeros(main_length, dtype=object)
i, j = n-1, m-1
xy = main_length-1
# ці змінні потрібні для того, щоб пам’ятати скільки було пропусків
gap_x = 0
gap_y = 0
while xy!=0:
    if trace[xy][0]=='M':
        seq_1_alignment[xy] = seq_1[i+gap_x]
        seq_2_alignment[xy] = seq_2[j+gap_y]
    elif trace[xy][0]=='By':
        seq_1_alignment[xy] = '-'
        seq_2_alignment[xy] = seq_2[j+gap_y]
        gap_x += 1
    elif trace[xy][0]=='Bx':
        seq_2_alignment[xy] = '-'
        seq_1_alignment[xy] = seq_1[i+gap_x]
        gap_y += 1
    i -= 1
    j -= 1
    xy -= 1

# створення гарних таблиць із результатами
multi_idx = pd.MultiIndex.from_product([seq_1, ['M', 'Bx', 'By']])
A = np.array([M, Bx, By])
F_df = pd.DataFrame(index=multi_idx, columns=seq_2)
for i in range(n):
    F_df.iloc[i*3,:] = M[i]
    F_df.iloc[i*3+1,:] = Bx[i]
    F_df.iloc[i*3+2,:] = By[i]
F_df = F_df.replace(-np.inf, "'-inf'")

T_df = pd.DataFrame(index=multi_idx, columns=seq_2)
for i in range(n):
    T_df.iloc[i*3,:] = M_T[i]
    T_df.iloc[i*3+1,:] = Bx_T[i]
    T_df.iloc[i*3+2,:] = By_T[i]

print('Completed matrix F:\n', F_df)
print('\nTraceback matrix T:\n', T_df)
print('\nAligned sequences:')
print(seq_1_alignment[1:])
print(seq_2_alignment[1:])

# збереження матриць до csv-файлів
F_df.to_csv('weights_affine.csv', sep=';')
T_df.to_csv('traceback_affine.csv', sep=';')