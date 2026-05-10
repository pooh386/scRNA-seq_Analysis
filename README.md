# scRNA-seq Analysis Pipeline: From Plant Science to Human Omics



## 1. Project Motivation
본 프로젝트는 식물생명공학 전공 배경의 연구자가 인체 오믹스(Human Omics) 데이터 분석으로 도메인을 확장하는 과정을 실증하기 위해 기획됨.

**도메인 확장성**: 식물 Wet-lab(mRNA 추출 및 Intron/Exon 분석) 경험을 바탕으로 데이터의 생물학적 본질을 이해하고 이를 데이터 분석 역량과 통합함.

**융합적 접근**: 위성 데이터(NDVI) 등 대규모 데이터 처리 경험을 단일 세포 전사체 데이터(scRNA-seq) 파이프라인 구축에 적용하여 Computational Biology 역량을 확보함



## 2. Detailed Analysis Workflow & Code


### Step 1: 데이터 품질 관리 (Quality Control)
분석의 정확도를 저해하는 저품질 세포 및 노이즈 유전자를 필터링함.
- 수행 내용: 3개 미만의 세포에서만 발현되는 유전자를 제거하여 시퀀싱 에러 및 통계적 유의성이 낮은 노이즈를 배제함.


```python
sc.pp.filter_genes(adata, min_cells=3)
```
<img width="1433" height="614" alt="스크린샷 2026-05-10 223435" src="https://github.com/user-attachments/assets/1dff9b6e-57de-4ffc-85fd-d594d08dddfa" />



### Step 2: 데이터 정규화 및 로그 변환 (Normalization)
세포별 시퀀싱 깊이(Depth) 차이를 보정하고 데이터 분포를 안정화함.
- 수행 내용: 모든 세포의 총 발현량을 동일 스케일($10^4$)로 맞춘 후, $log(x+1)$ 변환을 통해 데이터의 왜도(Skewness)를 줄임.

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```



### Step 3: 고변이 유전자 선택 (Feature Selection)
세포 군집화에 기여도가 높은 핵심 유전자를 식별함.
- 수행 내용: cell_ranger flavor를 사용하여 생물학적 변동성이 큰 유전자를 선별함. 이는 차원 축소의 효율성을 극대화하는 핵심 단계임.

```python
sc.pp.highly_variable_genes(adata, flavor='cell_ranger')
sc.pl.highly_variable_genes(adata)
```
<img width="1438" height="515" alt="스크린샷 2026-05-10 223714" src="https://github.com/user-attachments/assets/4993929f-8d40-4be6-8a43-ddefa0cd3d16" />



### Step 4: 차원 축소 및 이웃 그래프 생성 (PCA & Neighbors)
고차원의 데이터를 정보 손실을 최소화하며 저차원으로 압축함.
-수행 내용: PCA를 통해 주성분을 추출하고, 이를 기반으로 세포 간의 거리(유사도)를 계산하여 KNN(K-Nearest Neighbors) 그래프를 구축함.

```python
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
```



### Step 5: UMAP 시각화 및 군집화 (Clustering)
세포들 간의 관계를 2차원 평면에 시각화하고 통계적으로 군집화함.
- 수행 내용: UMAP 알고리즘을 통한 좌표 임베딩을 수행하고, Leiden 알고리즘을 적용하여 세포들을 유의미한 군집으로 분류함.

```python
sc.tl.umap(adata)
sc.tl.leiden(adata, flavor='igraph', n_iterations=2)
sc.pl.umap(adata, color='leiden')
```
<img width="824" height="612" alt="스크린샷 2026-05-10 223651" src="https://github.com/user-attachments/assets/9f359998-f757-4dac-bfd4-51395271a665" />



### Step 6: 마커 유전자 분석 (Cell Type Annotation)
각 군집의 생물학적 정체성을 검증함.
- 수행 내용: 주요 면역 세포 마커(CST3, NKG7, MS4A1) 발현량을 UMAP에 매핑하여 B세포, T세포 등의 분포를 식별함.

```python
sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7', 'MS4A1'])
```
<img width="1444" height="400" alt="스크린샷 2026-05-10 223732" src="https://github.com/user-attachments/assets/f1d9e9d5-9888-411d-ad9f-b455ba218e39" />




## 3. Troubleshooting & Engineering Experience
환경 구축 과정에서 발생한 기술적 이슈를 주도적으로 해결함.

1. Python 버전 최적화: 최신 Python 버전에서의 라이브러리 의존성(numba, leidenalg) 충돌을 확인하고, Python 3.11 기반 가상환경 구축을 통해 분석 안정성을 확보함.
2. 저장소 경량화: 100MB 초과 이진 파일(llvmlite.dll)에 의한 Push 에러 발생 시, .gitignore 설정 및 Git 캐시 리셋을 수행하여 분석 코드 중심의 클린한 저장소를 유지함.



## 4. How to Use
1. Clone this repository.
2. Build environment: py -3.11 -m venv venv
3. Install dependencies: pip install scanpy leidenalg igraph matplotlib
4. Run Jupyter Notebook: PBMC3k_Analysis.ipynb

