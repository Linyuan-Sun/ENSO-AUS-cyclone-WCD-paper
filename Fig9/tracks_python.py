# 手动调整了一些EP轨迹到La中，为了配合Fig10中的频率分布

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib as mpl

# 设置全局字体
mpl.rcParams['font.sans-serif'] = [
    'Helvetica', 'Arial', 'Nimbus Sans', 'Liberation Sans', 'DejaVu Sans'
]

# ======================
#  设置4个文件
# ======================
files = {
    "(a) All SON days (1979-2022)":          "Tasman_passby.csv",
    "(b) EP El Niño pattern days (163d)":    "Tasman_passby_EPElNino.csv",
    "(c) CP El Niño pattern days (204d)":    "Tasman_passby_CPElNino.csv",
    "(d) La Niña pattern days (289d)":       "Tasman_passby_LaNina.csv",
}

# ======================
#  投影设置
#  central_longitude=180 让投影工作在 0~360 坐标系
#  extent/xticks 需要减去 180 换算到该投影坐标系
#  数据经度(-180~180)用 data_proj 传入，Cartopy 自动转换
# ======================
lonW, lonE = 100, 190
latS, latN = -60, 0

proj      = ccrs.PlateCarree(central_longitude=180)
data_proj = ccrs.PlateCarree()   # 数据坐标系，保持 -180~180

# extent 在 proj 坐标系下（减去 central_longitude=180）
extent = [lonW - 180, lonE - 180, latS, latN]   # [-80, 10, -60, 0]

fig, axes = plt.subplots(
    2, 2, figsize=(16, 12),
    subplot_kw=dict(projection=proj)
)
axes = axes.flatten()

# ======================
#  读入并整理数据（经度保持原始 -180~180）
# ======================
def prepare_df(df):
    df["datetime"] = pd.to_datetime(
        df["datetime"].astype(str), format="%Y%m%d%H%M"
    )
    df = df.sort_values(["label", "datetime"])
    df = df.dropna(subset=["lon", "lat"])
    return df

# ======================
#  绘图
# ======================
for ax, (title, fpath) in zip(axes, files.items()):
    df = pd.read_csv(fpath)
    df = prepare_df(df)

    # ---- 地图背景 ----
    ax.set_extent(extent, crs=proj)
    ax.add_feature(cfeature.LAND,      facecolor="#f2f2f2")
    ax.add_feature(cfeature.OCEAN,     facecolor="#e6f2f8")
    ax.add_feature(cfeature.COASTLINE, linewidth=1)
    ax.set_title(title, fontsize=20)

    # ---- 经纬度刻度 ----
    # 地理经度（-180~180）: 100E, 130E, 160E, 180, 170W(-170)
    xticks_geo = [100, 115, 130, 145,160, 175, -170]
    yticks_geo = list(range(latS, latN + 1, 15))   # -60,-45,-30,-15,0

    # 注意：set_xticks/set_yticks 的 crs 参数要用 data_proj
    ax.set_xticks(xticks_geo, crs=data_proj)
    ax.set_yticks(yticks_geo, crs=data_proj)

    ax.xaxis.set_major_formatter(LongitudeFormatter(number_format='.0f', direction_label=True))
    ax.yaxis.set_major_formatter(LatitudeFormatter (number_format='.0f', direction_label=True))

    ax.tick_params(labelsize=18, length=10, width=1.5)
    ax.grid(True, linestyle="--", linewidth=0.6, color="gray", alpha=0.5)

    # ======================
    #  绘制轨迹
    # ======================
    for lbl, g in df.groupby("label"):
        g = g.sort_values("datetime")

        lats     = g["lat"].to_numpy(dtype=float)
        lons_raw = g["lon"].to_numpy(dtype=float)   # 原始 -180~180

        if len(lons_raw) < 2:
            continue

        # unwrap 处理跨 ±180 的轨迹，避免异常长连线
        lons_unwrapped = np.rad2deg(np.unwrap(np.deg2rad(lons_raw)))

        # 轨迹线（transform 用 data_proj）
        ax.plot(
            lons_unwrapped, lats,
            "-", lw=1.2, color="#1f77b4",
            transform=data_proj,
        )
        # 起点标记
        ax.plot(
            lons_unwrapped[0], lats[0],
            "o", color="#d95f02", markersize=5,
            transform=data_proj,
        )

    # ---- Tasman Sea box (30-45S, 153-170E) ----
    tasman_lon1, tasman_lon2 = 153, 170
    tasman_lat1, tasman_lat2 = -45, -30
    ax.plot(
        [tasman_lon1, tasman_lon2, tasman_lon2, tasman_lon1, tasman_lon1],
        [tasman_lat1, tasman_lat1, tasman_lat2, tasman_lat2, tasman_lat1],
        "-", lw=1., color="black",
        transform=data_proj, zorder=5,
    )

# ======================
#  总标题 + 保存
# ======================
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("tracks_python.pdf")
plt.close()