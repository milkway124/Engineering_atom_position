## Engineering_atom_position

Density Functional Theory이 적용된 [Quantum espresso](https://www.quantum-espresso.org/)(QE)와 [VASP](https://www.vasp.at/wiki/index.php/The_VASP_Manual)에 계산하기 위해서는 atom들이 위치가 필수적이다.

이 atom들의 위치는 POSCAR 파일 (Atom의 위치와 lattice 구조가 담긴 파일)에 존재하며, 이 atom들은
fractional coordinate로 기술 되어있다. 이 atom들을 원하는 방향으로 변형시키는 코드들을 총 5가지에 대해 알아본다.

관련된 모든 코드는 linux에서 사용을 하기 위해서 작성이 되었습니다.

또한 언어는 python과 bash언어로 이루어져 있습니다.
___
#### 1. 작은 unit-cell 생성 및 예제
기본적인 구조 파일 (POSCAR 파일)은 [VESTA](https://jp-minerals.org/vesta/en/)을 이용하여 시각적으로 파악한다.

Example : Diamond, Silicon, h-BN ...


___
#### 2. 작은 unit-cell을 큰 크기의 supercell로 늘리기
작은 unit-cell을 [VESTA](https://jp-minerals.org/vesta/en/)을 이용하여 늘릴 수도 있지만,
직접 작성한 코드를 이용하여 원하는 크기로 늘리는 작업을 수행한다.


___
#### 3. Defect를 supercell에 넣기
크게 키운 supercell에 원하는 defect을 넣는다.

Example : VB 결함이 있는 h-BN, VNNB 결함이 있는 Diamond, ...

추가적으로 Defect 주위에 atom들도 다른 이름으로 특정을 짓어 POSCAR 파일을 재생성한다.


___
#### 4. 2차원 supercell을 구조 변형 (Bubble, Cone, Pillar, Cross 모양)
2차원 물질중에 대표적으로 h-BN을 이용하여 구조 변형을 시킨다.


___
#### 5. 2차원 supercell을 Stacking하기
2차원 물질은 단층으로 이루어진 물질이 아닌 3차원으로 쌓여있기 때문에, 
코드를 이용하여 복층을 구성한다.


___
#### 번외1. POSCAR 파일에 있는 atom의 position을 fix하기 (0 0 0)



#### 번외2. POSCAR 파일을 QE의 input-file (scf.in)으로 넣기



#### 번외3. QE을 통해 Relaxed된 파일을 POSCAR 파일로 변환하기


