## Engineering_atom_position

Density Functional Theory이 적용된 Quantum espresso(QE)와 VASP에 사용하기 위해,
POSCAR 파일 (Atom의 위치와 lattice 구조가 담긴 파일)을 원하는 방향으로 변형시키는 코드들

### 1. 작은 unit-cell
기본적인 구조 파일 (POSCAR 파일)은 [VESTA](https://jp-minerals.org/vesta/en/)을 이용하여 눈으로 파악한다.

### 2. supercell로 늘리기

### 3. Defect를 supercell에 넣기

### 4. 2차원 supercell을 Stacking하기

### 번외1. POSCAR 파일에 있는 atom의 position을 fix하기 (0 0 0)

### 번외2. POSCAR 파일을 QE의 input-file (scf.in)으로 넣기

### 번외3. QE을 통해 Relaxed된 파일을 POSCAR 파일로 변환하기
