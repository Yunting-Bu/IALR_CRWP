# 动能算符对比：schprog1.f（两区域） vs IALR_CRWP（三区域）

## schprog1.f 的做法

每个角量子数通道 `I=1,NAP`，在一个大的 `!$OMP DO` 循环中完成：

| 步骤 | 操作 | 对象 | DST大小 |
|---|---|---|---|
| 1 | PLANI(1,I) + CENT3 + PLANI(2,I) | Z动能，r > NRASY 的额外基 | DST(**NZINT**) x (NR-NRASY)个 |
| 2 | PLANA(1,I) + CENT1 + PLANA(2,I) | Z动能，r <= NRASY 的共享基 | DST(**NZ**) x NRASY个 |
| 3 | PLANI(3,I) + CENT2 + PLANI(4,I) | r动能，相互作用区 | DST(**NR**) x NZINT个 |
| 4 | DGEMM + CENT4 + DGEMM | r动能，渐近区 | DGEMM(NZ x NRASY^2) |

**关键特点**：共享的 NRASY 个 r 基函数，Z 动能在**整个 NZ 网格**上做一次长 DST，结果按 Z 范围拆分到两个区域。

## IALR_CRWP（三区域）的做法

`kinAction()` 分步执行：

| 步骤 | 操作 | 对象 | DST大小 |
|---|---|---|---|
| 1 | copy int->asy->lr | 把波包数据复制到各区域 | -- |
| 2 | dsint + kinEigen + dsint | 相互作用区 Z+r 动能一起做 | DST(**nZ_I=62**) x vint个 |
| 3 | dsint + kinEigen + dsint | 渐近区 Z+r 动能一起做 | DST(**nZ_IA=295**) x vasy个 |
| 4 | dsint + kinEigen + dsint | 长程区 Z+r 动能一起做 | DST(**nZ_IALR=2047**) x vlr个 |
| 5 | copy lr->asy->int | 结果复制回来 | -- |
| 6 | daxpy 累加 | 加到 TDWP 上 | -- |

**关键特点**：`kinEigen(iZ,ir) = ZkinMat(iZ) + rKinMat(ir)`，因为 r 用 PO-FBR 表示，r 动能已经是对角的，可以**和 Z 动能合并在一次 DST 中完成**，不需要单独的 r 动能步骤。

## 计算量对比（FH2 参数）

以每个角通道为单位（NZ=2047, NZINT=62, NR=100, NRASY=25）：

**schprog1.f 每个角通道**：
```
Z动能:  2 x 25 x DST(2047) + 2 x 75 x DST(62)   <-- 25个长DST是主要开销
r动能:  2 x 62 x DST(100) + 2 x DGEMM(~1985 x 25 x 25)
```

**IALR_CRWP 每个等效角通道**：
```
Z+r动能: 2 x 100 x DST(62) + 2 x 25 x DST(295) + 2 x 1 x DST(2047)
复制:    O(nZ_I x vasy) 两次
```

把 DST 代价按 `N log N` 估算：

| | schprog1 | IALR_CRWP |
|---|---|---|
| 长 DST | 25 x DST(2047) ~ 25 x 22,500 = **562,500** | 1 x DST(2047) ~ **22,500** |
| 中 DST | -- | 25 x DST(295) ~ 25 x 2,400 = **60,000** |
| 短 DST | 75 x DST(62) ~ 75 x 370 = **27,750** | 100 x DST(62) ~ **37,000** |
| r 动能 | 62 x DST(100) + DGEMM ~ **290,000** | **0**（已合并） |
| **合计** | **~880,000** | **~120,000** |

## 结论

**IALR_CRWP 三区域方案在动能算符上效率大约高 7 倍**，原因有两个：

1. **避免了长 DST 的重复计算**。schprog1 对所有 NRASY=25 个共享 r 基函数都做 DST(2047)，而 IALR_CRWP 只对 vlr=1 个 r 基函数做 DST(2047)，其余 25 个只做 DST(295)。这是最大的节省——DST(2047) 比 DST(295) 贵约 **9 倍**，而做的次数从 25 降到了 1。

2. **r 动能零额外开销**。IALR_CRWP 用 PO-FBR 表示 r，动能已经是对角的，`kinEigen(iZ,ir) = TZ(iZ) + Tr(ir)` 在 Z 动量空间中一步完成。schprog1 用 particle-in-box DVR，需要额外的 r 方向 DST 和渐近区 DGEMM。

**schprog1 唯一的优势**：它用 FFTW（`DFFTW_EXECUTE`），比 DFFTPACK（`dsint`）在同等大小下快 2-3 倍左右。但这不足以抵消算法上 7 倍的优势。如果把 `dsint` 换成 FFTW，效率优势会更大。

## 如果 schprog1 也用三区域？

如果 schprog1 也拆成三个区域（相互作用、渐近、长程），Z 动能的 DST 开销将**完全一致**。剩下的差别只在 **r 动能**的处理方式。

### Z 动能（两者相同）

| 区域 | DST 开销 |
|---|---|
| 相互作用 | 2 x 100 x DST(62) ~ 37,000 |
| 渐近 | 2 x 25 x DST(295) ~ 60,000 |
| 长程 | 2 x 1 x DST(2047) ~ 22,500 |
| **Z 动能合计** | **~119,500** |

### r 动能（只有 schprog1 需要）

| 区域 | 操作 | 开销 |
|---|---|---|
| 相互作用 | 2 x nZ_I x DST(NR) | 2 x 62 x DST(100) ~ 86,800 |
| 渐近 | 2 x DGEMM(233 x 25 x 25) | ~ 291,000 |
| 长程 | 忽略 (vlr=1) | ~0 |
| **r 动能合计** | | **~378,000** |

### 对比

| | schprog1 (三区域假设) | IALR_CRWP |
|---|---|---|
| Z 动能 | 119,500 | 119,500 |
| r 动能 | **378,000** | **0** |
| **合计** | **~497,500** | **~119,500** |

优势从 ~7 倍缩小到 **~4 倍**，完全来自 **PO-FBR 表示下 r 动能免费**这一点。

### 根本原因：r 基函数表示不同

| | schprog1 | IALR_CRWP |
|---|---|---|
| r 表示 | particle-in-box DVR (正弦基) | PO-FBR (势优化基) |
| r 动能 | 需要 DST 变换到动量空间 | **已经对角**，eigenvalue = PO 本征值 |
| r 势能 | 已经在 DVR，直接乘 | 需要 rBasisTrans (FBR<->DVR) |

IALR_CRWP 把 r 动能的代价**转移到了势能步骤**（通过 `rBasisTrans` 做 FBR<->DVR 变换）。但这并不亏，因为：

- 势能步骤中的 `rBasisTrans` 无论如何都要做角度方向的 DGEMM，r 变换只是多乘一次矩阵
- 动能步骤在时间步循环中调用频率和势能一样（每步都调），省掉 r 动能的 DST/DGEMM 是**净节省**

因此即使都用三区域，**PO-FBR 方案仍然有 ~4 倍优势**。
