# ParityGate

一个用于量子计算的Julia包，实现了奇偶门（Parity Gate）功能。

## 安装

```julia
using Pkg
Pkg.add("https://github.com/yourusername/parity_gate.git")
```

## 使用方法

```julia
using ParityGate
using Yao

# 创建一个2量子比特的奇偶门
circuit = parity_gate(2)

# 创建一个量子寄存器
reg = zero_state(2)

# 应用奇偶门
apply_parity_gate!(reg, circuit)
```

## 功能

- `parity_gate(n::Int)`: 创建一个n量子比特的奇偶门
- `apply_parity_gate!(reg::ArrayReg, circuit::ChainBlock)`: 将奇偶门应用到量子寄存器上

## 示例

```julia
# 创建一个3量子比特的奇偶门
circuit = parity_gate(3)

# 创建一个初始状态|111⟩
reg = product_state(3, 0b111)

# 应用奇偶门
apply_parity_gate!(reg, circuit)
```

## 许可证

MIT 