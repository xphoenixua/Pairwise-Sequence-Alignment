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
    F[i, 0] = -i * d
    T[i, 0] = tuple((i, 0))
for j in range(m):
    F[0, j] = -j * d
    T[0, j] = tuple((0, j))

# заповнення комірок матриць
for i in range(1, n):
    for j in range(1, m):
        xy = str(seq_1[i]+seq_2[j])
        max_arr = [F[i-1, j-1]+s[xy], F[i-1, j]-d, F[i, j-1]-d]
        F[i, j] = np.max(max_arr)
        # залежно від варіанту запам’ятовую звідки було здійснено переміщення
        idx = np.argmax(max_arr)
        if idx==0:
            T[i, j] = tuple((i-1, j-1))
        elif idx==1:
            T[i, j] = tuple((i-1, j))
        elif idx==2:
            T[i, j] = tuple((i, j-1))

# відстеження зворотнього шляху від комірки F(n,m) до комірки F(0,0)
trace = np.zeros(m+1, dtype=tuple)
# наступний рядок потрібен як початкова умова для подальшого вирівнювання
trace[0] = (6,6)
i, j = n-1, m-1
xy = 1
while i!=0 and j!=0:
    trace[xy] = T[i, j]
    i, j = T[i, j][0], T[i, j][1]
    xy += 1

# вирівнювання двох послідовностей
seq_1_alignment = np.zeros(n, dtype=object)
seq_2_alignment = np.zeros(m, dtype=object)
i, j = n-1, m-1
# ці змінні потрібні для того, щоб пам’ятати скільки було пропусків
gap_x = 0
gap_y = 0
for xy in range(1, trace.shape[0]):
    if np.abs(trace[xy][0]-trace[xy-1][0])>0:
        seq_1_alignment[i] = seq_1[i+gap_x]
        i -= 1
    elif np.abs(trace[xy][0]-trace[xy-1][0])==0:
        seq_1_alignment[i] = '-'
        gap_x += 1
        i -= 1
    if np.abs(trace[xy][1]-trace[xy-1][1])>0:
        seq_2_alignment[j] = seq_2[j+gap_y]
        j -= 1
    elif np.abs(trace[xy][1]-trace[xy-1][1])==0:
        seq_2_alignment[trace[xy][1]] = '-'
        gap_y += 1
        j -= 1

# виведення остаточних результатів
score = F[n-1, m-1]
F_df = pd.DataFrame(data=F, index=seq_1, columns=seq_2)
print('Completed matrix F:\n', F_df)
T_df = pd.DataFrame(data=T, index=seq_1, columns=seq_2)
print('\nTraceback matrix T:\n', T_df)
print('\nScore: ', score)
print('Aligned sequences:')
print(seq_1_alignment)
print(seq_2_alignment)

# збереження матриць до csv-файлів
F_df.to_csv('F_global.csv', sep=';')
T_df.to_csv('T_global.csv', sep=';')