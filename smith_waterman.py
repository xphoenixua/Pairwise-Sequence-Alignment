import numpy as np
import pandas as pd
# даний модуль уже містить матриці а-к замін
import blosum as bl

# початкові умови
d = 8
s = bl.BLOSUM(62)
seq_1 = 'FLEKDG'
seq_2 = 'IQLEKK'
seq_1 = [0] + [*seq_1]
seq_2 = [0] + [*seq_2]

# ініціалізація матриці ваг та матриці зворотнього проходу
n, m = len(seq_1), len(seq_2)
F = np.zeros((n, m))
T = np.zeros((n, m), dtype=tuple)
for i in range(n):
    F[i, 0] = 0
    T[i, 0] = tuple((0,0))
for j in range(m):
    F[0, j] = 0
    T[0, j] = tuple((0,0))

# заповнення комірок матриць
for i in range(1, n):
    for j in range(1, m):
        xy = str(seq_1[i]+seq_2[j])
        max_arr = [0, F[i-1, j-1]+s[xy], F[i-1, j]-d, F[i, j-1]-d]
        F[i, j] = np.max(max_arr)
        # залежно від варіанту запам’ятовую звідки було здійснено переміщення
        idx = np.argmax(max_arr)
        if idx==0:
            T[i, j] = tuple((0, 0))
        elif idx==1:
            T[i, j] = tuple((i-1, j-1))
        elif idx==2:
            T[i, j] = tuple((i-1, j))
        elif idx==3:
            T[i, j] = tuple((i, j-1))

# визначаємо максимальне значення зі всіх комірок матриці
# і запам’ятовуємо її координати
score = np.max(F)
score_idx = np.argmax(F)
r_idx = score_idx // F.shape[1]
c_idx = score_idx % F.shape[1]

# відстеження зворотнього шляху від комірки max(F(i,j)) до найближчої комірки
# комірки зі значенням (0,0)
trace = np.zeros(c_idx, dtype=tuple)
i, j = r_idx, c_idx
xy = 0
while i!=0 and j!=0:
    trace[xy] = T[i, j]
    i, j = T[i, j][0], T[i, j][1]
    xy += 1

# вирівнювання двох послідовностей
seq_1_alignment = np.zeros(r_idx, dtype=object)
seq_2_alignment = np.zeros(c_idx, dtype=object)
i, j = r_idx, c_idx
# ці змінні потрібні для того, щоб пам’ятати скільки було пропусків
gap_x = 0
gap_y = 0
for xy in range(1, c_idx):
    if trace[xy]==0:
        break
    if np.abs(trace[xy][0]-trace[xy-1][0])>0:
        seq_1_alignment[i-1] = seq_1[i+gap_x]
        i -= 1
    elif np.abs(trace[xy][0]-trace[xy-1][0])==0:
        seq_1_alignment[i-1] = '-'
        gap_x += 1
        i -= 1
    if np.abs(trace[xy][1]-trace[xy-1][1])>0:
        seq_2_alignment[j-1] = seq_2[j+gap_y]
        j -= 1
    elif np.abs(trace[xy][1]-trace[xy-1][1])==0:
        seq_2_alignment[j-1] = '-'
        gap_y += 1
        j -= 1

seq_1_alignment = seq_1_alignment[seq_1_alignment!=0]
seq_2_alignment = seq_2_alignment[seq_2_alignment!=0]

# виведення остаточних результатів
F_df = pd.DataFrame(data=F, index=seq_1, columns=seq_2)
print('Completed matrix F:\n', F_df)
T_df = pd.DataFrame(data=T, index=seq_1, columns=seq_2)
print('\nTraceback matrix T:\n', T_df)
print('\nScore: ', score)
print('Aligned sequences:')
print(seq_1_alignment)
print(seq_2_alignment)

# збереження матриць до csv-файлів
F_df.to_csv('F_local.csv', sep=';')
T_df.to_csv('T_local.csv', sep=';')