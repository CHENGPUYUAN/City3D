# CLAUDE.md — City3D 开发路由

City3D:机载 LiDAR 点云 → 建筑三维模型重建(假设生成 + 选择优化)。
基于 PolyFit,核心优化是 binary linear program 的求解。

## 快速命令

```bash
# 构建(本地已配置好 build/)
cmake --build build --target CLI_Example_1 -j8

# 重建(GUI 可执行文件,需 Qt)
cmake --build build --target City3D -j8

# 运行(CLI;第 6 个参数可选求解器: cuopt|gurobi|mindoptpy|highs|scip,缺省 auto=cuopt 优先、失败回落本地链)
./build/bin/CLI_Example_1 <cloud.ply|las|laz> <footprint.obj> <out.obj> [pixel_size] [min_points] [solver]

# 测试数据
./build/bin/CLI_Example_1 data/002.ply data/002_footprint.obj /tmp/out.obj
```

Python 侧(数据转换)用 `uv run`,见 [[build-and-run]]。

## 目录速览

| 路径 | 职责 |
|---|---|
| `code/method/` | 重建管线:分割 → 屋面平面提取 → 假设生成 → face selection 优化 |
| `code/math/` | LinearProgram 抽象 + 各求解器后端(CUOPT / GUROBI / HIGHS / MINDOPTPY / SCIP) |
| `code/model/` | PointSet / Map 数据结构与 IO |
| `code/CLI_Example_1/` | 单场景 CLI 入口 |
| `code/3rd_party/3rd_scip/` | bundled SCIP 8.0.3 备用(系统装了 scipoptsuite 10 时优先动态链接系统库) |
| `data/dsm_to_city3d.py` | DSM 栅格 + footprint GeoJSON → City3D 输入(PLY/OBJ) |

## 文档索引

按 scope/tags 路由,先读 frontmatter 再决定是否读全文:

- **核心架构** `docs/core/`
  - [[solver-backends]] — LP 求解器分发:GUROBI → HIGHS → MINDOPTPY → SCIP 策略、双 SCIP 符号冲突、MindOpt 社区许可证
  - [[reconstruction-pipeline]] — 重建管线各阶段、中间产物、失败路径
- **决策记录** `docs/decisions/`
  - [[use-mindoptpy-backend]] — 为什么用 mindoptpy 子进程取代 gurobipy(社区许可证无规模限制)
  - [[reject-parallel-scip-ug]] — 为什么不用并发 SCIP / FiberSCIP / UG(含实测数据)
- **运维** `docs/ops/`
  - [[build-and-run]] — 构建细节(许可门控、uv 依赖)、DSM 数据管线、常见坑

## 关键约束(改动前必读)

1. **MindOpt 许可证鉴权**:MINDOPTPY 后端走阿里云 license server(需联网,
   ~0.2s/次);SDK 两处 `libmindopt.so` 的 RWE GNU_STACK 已打补丁,重装 MindOpt
   需重新清 PT_GNU_STACK 的 PF_X 位(见 [[solver-backends]])。
2. **Logger 陷阱**:CLI 里调 `Logger::instance()->set_value(...)` 前必须先
   `Logger::initialize()`,否则空指针段错误(GUI 不受影响,它在 main_window 里先初始化了)。
3. **候选面上限**:`Method::max_allowed_candidate_faces = 30000`(method_global.cpp),
   超过的建筑直接跳过。DSM 高密度数据常触顶(见 [[build-and-run]] 的 case1 说明)。
4. `code/3rd_party/` 是裁剪过的第三方树,不要整体升级;对它的改动(如 TPI 启用)都有本地原因。

上游通用文档:ReadMe.md(项目介绍、依赖)、Parameters.md(参数调优)。
