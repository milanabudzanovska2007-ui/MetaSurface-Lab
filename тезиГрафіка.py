"""MetaSurface Lab — скалярний симулятор із порівнянням фазового квантування.

Реалізує запропонований ІБВП:
D_B = 1 - (1 - E_i) * (1 - E_phi) * (1 - E_eta)

E_i   — нормована L1-похибка інтенсивності;
E_phi — зважена фазова похибка після вилучення глобальної фази;
E_eta — абсолютна різниця часток інтенсивності в цільовій області.

Модель навчальна: це скалярне поширення методом кутового спектра,
а не повновекторне розв'язання рівнянь Максвелла чи модель наногеометрії.
"""

from __future__ import annotations

import csv
import math
import argparse
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
from scipy.fft import fft2, ifft2, fftfreq
from PIL import Image, ImageTk

APP_NAME = "MetaSurface Lab — ІБВП"


def wrap_phase(phi: np.ndarray) -> np.ndarray:
    """Загортає фазу в інтервал [-pi, pi)."""
    return np.angle(np.exp(1j * phi))


def make_phase_profile(kind: str, x: np.ndarray, y: np.ndarray,
                       wavelength: float, focal_length: float,
                       axicon_angle_deg: float, vortex_charge: int,
                       random_strength: float, seed: int) -> np.ndarray:
    """Створює фазовий профіль у радіанах; фізичні координати — метри."""
    k = 2.0 * np.pi / wavelength
    r = np.sqrt(x * x + y * y)

        # Лінза змінює фазу залежно від відстані до центра.
    if kind == "Лінза":
        phi = -k * (x * x + y * y) / (2.0 * focal_length)
        # Аксикон задає фазу, що залежить від радіальної відстані до центра.
    elif kind == "Аксикон":
        theta = np.deg2rad(axicon_angle_deg)
        phi = -k * r * np.sin(theta)
        # Вихрова пластинка створює фазовий вихор; m задає його заряд.
    elif kind == "Вихрова пластинка":
        phi = int(vortex_charge) * np.arctan2(y, x)
        # Лінійна зміна фази моделює нахил хвильового фронту.
    elif kind == "Керований нахил":
        period = max(8.0 * wavelength, min(abs(focal_length), 1e-3))
        phi = 2.0 * np.pi * x / period
        # Тут фаза випадкова, але seed робить результат відтворюваним.
    elif kind == "Випадкова фаза":
        rng = np.random.default_rng(int(seed))
        phi = rng.uniform(-np.pi, np.pi, size=x.shape) * float(random_strength)
    else:
        phi = np.zeros_like(x)
    return wrap_phase(phi)


def quantize_phase(phase: np.ndarray, levels: int) -> np.ndarray:
    """Квантує фазу до рівномірних рівнів на колі; levels має бути >= 2."""
    if int(levels) < 2:
        raise ValueError("Кількість фазових рівнів має бути не меншою за 2.")
    step = 2.0 * np.pi / int(levels)
    # Найближчий рівень на періодичному інтервалі [-pi, pi).
    indices = np.mod(np.floor((wrap_phase(phase) + np.pi) / step + 0.5), int(levels))
    qphase = -np.pi + indices * step
    return wrap_phase(qphase)


def angular_spectrum(field: np.ndarray, wavelength: float, dx: float,
                     z: float) -> np.ndarray:
    """Поширення комплексного поля методом кутового спектра."""
    ny, nx = field.shape
    fx = fftfreq(nx, d=dx)
    fy = fftfreq(ny, d=dx)
    FX, FY = np.meshgrid(fx, fy)
        # q визначає, чи є просторова частота поширюваною в межах цієї скалярної моделі.
    q = 1.0 / wavelength**2 - FX**2 - FY**2
    transfer = np.zeros(field.shape, dtype=np.complex128)
        # Залишаємо лише компоненти, для яких корінь у передавальній функції дійсний.
    propagating = q >= 0
    transfer[propagating] = np.exp(1j * 2.0 * np.pi * z * np.sqrt(q[propagating]))
        # Перетворюємо поле у частотну область, застосовуємо поширення і повертаємося назад.
    return ifft2(fft2(field) * transfer)


def _energy_fraction(intensity: np.ndarray, roi_fraction: float = 0.20) -> float:
    """Частка суми інтенсивності в центральній квадратній ROI."""
    n_y, n_x = intensity.shape
    half_y = max(1, int(round(n_y * roi_fraction / 2)))
    half_x = max(1, int(round(n_x * roi_fraction / 2)))
    cy, cx = n_y // 2, n_x // 2
    roi = intensity[max(0, cy-half_y):min(n_y, cy+half_y),
                    max(0, cx-half_x):min(n_x, cx+half_x)]
    total = float(np.sum(intensity))
    return float(np.sum(roi) / total) if total > 0 else 0.0


def calculate_index(field_ref: np.ndarray, field_test: np.ndarray,
                    roi_fraction: float = 0.20) -> dict[str, float]:
    """Обчислює три складові й ІБВП для полів однакової форми.

    Якщо обидва поля нульові, E_phi=0. Якщо нульове лише одне поле,
    E_phi=1. Це явна домовленість для невизначеної фази.
    """
    if field_ref.shape != field_test.shape:
        raise ValueError("Еталонне й тестове поля мають мати однакову форму.")
    if not (0.0 < roi_fraction <= 1.0):
        raise ValueError("Частка ROI має бути в інтервалі (0, 1].")

        # Інтенсивність — квадрат модуля комплексного поля; саме її показують карти яскравості.
    i_ref = np.abs(field_ref) ** 2
    i_test = np.abs(field_test) ** 2
    denom_i = float(i_ref.sum() + i_test.sum())
    if denom_i == 0.0:
        e_i = 0.0
    else:
                # E_i порівнює інтенсивності в усіх пікселях і нормує різницю на сумарну інтенсивність.
        e_i = float(np.sum(np.abs(i_ref - i_test)) / denom_i)
    e_i = float(np.clip(e_i, 0.0, 1.0))

    norm_ref = float(np.sum(np.abs(field_ref) ** 2))
    norm_test = float(np.sum(np.abs(field_test) ** 2))
    if norm_ref == 0.0 and norm_test == 0.0:
        e_phi = 0.0
    elif norm_ref == 0.0 or norm_test == 0.0:
        e_phi = 1.0
    else:
        # theta прибирає глобальний фазовий зсув тестового поля відносно еталона.
                # Перекриття полів використовується, щоб оцінити спільний глобальний фазовий зсув.
        overlap = np.sum(np.conj(field_ref) * field_test)
        theta = float(np.angle(overlap))
                # Віднімаємо глобальний зсув: нас цікавить локальна різниця фаз, а не спільний зсув усього поля.
        phase_difference = np.angle(field_test * np.conj(field_ref) * np.exp(-1j * theta))
                # Точки, де амплітуда одного з полів мала, отримують меншу вагу у фазовому порівнянні.
        weights = np.abs(field_ref) * np.abs(field_test)
        weight_sum = float(weights.sum())
        e_phi = (float(np.sum(weights * np.abs(phase_difference)) / (np.pi * weight_sum))
                 if weight_sum > 0 else 1.0)
    e_phi = float(np.clip(e_phi, 0.0, 1.0))

        # E_eta порівнює частки інтенсивності, що залишилися в центральній області ROI.
    eta_ref = _energy_fraction(i_ref, roi_fraction)
    eta_test = _energy_fraction(i_test, roi_fraction)
    e_eta = float(np.clip(abs(eta_ref - eta_test), 0.0, 1.0))

        # ІБВП об’єднує три складові в одне число: більше значення означає більше сумарне відхилення.
    d_b = 1.0 - (1.0 - e_i) * (1.0 - e_phi) * (1.0 - e_eta)
    d_b = float(np.clip(d_b, 0.0, 1.0))
    return {"E_i": e_i, "E_phi": e_phi, "E_eta": e_eta,
            "D_B": d_b, "eta_ref": eta_ref, "eta_test": eta_test}


def simulate_pair(kind: str = "Лінза", wavelength_nm: float = 532.0,
                  pixel_um: float = 4.0, size: int = 256,
                  focal_mm: float = 20.0, z_mm: float = 20.0,
                  beam_radius_um: float = 150.0, levels: int = 8,
                  axicon_angle_deg: float = 1.0, vortex_charge: int = 1,
                  random_strength: float = 1.0, seed: int = 2026,
                  aperture_enabled: bool = True,
                  roi_fraction: float = 0.20) -> dict:
    """Рахує еталонне поле з неперервною фазою й тестове після квантування."""
    if wavelength_nm <= 0 or pixel_um <= 0 or focal_mm <= 0 or z_mm < 0 or beam_radius_um <= 0:
        raise ValueError("λ, крок пікселя, фокус і радіус мають бути додатними; z не може бути від'ємним.")
    if size not in (128, 256, 512, 1024):
        raise ValueError("Розмір сітки має бути 128, 256, 512 або 1024.")
    if levels < 2:
        raise ValueError("Кількість фазових рівнів має бути не меншою за 2.")
    if not (0 <= axicon_angle_deg <= 30) or random_strength < 0:
        raise ValueError("Кут аксикона має бути 0–30°, множник випадкової фази — невід'ємним.")

    wavelength = wavelength_nm * 1e-9
    dx = pixel_um * 1e-6
    focal_length = focal_mm * 1e-3
    z = z_mm * 1e-3
        # Створюємо координатну сітку в метрах; dx переводить крок із мікрометрів у SI.
    coords = (np.arange(size) - size / 2) * dx
    x, y = np.meshgrid(coords, coords)
    r = np.sqrt(x*x + y*y)
    phase_ref = make_phase_profile(kind, x, y, wavelength, focal_length,
                                   axicon_angle_deg, vortex_charge,
                                   random_strength, seed)
    w0 = max(beam_radius_um * 1e-6, dx)
        # Гаусів пучок найінтенсивніший біля центра і плавно слабшає до країв.
    amplitude = np.exp(-(r*r) / (w0*w0))
    if aperture_enabled:
        amplitude *= (r <= size * dx * 0.47)
    input_field = amplitude.astype(np.complex128)

        # Еталон використовує неперервну фазу, а тестовий варіант — її дискретизовану версію.
    phase_test = quantize_phase(phase_ref, levels)
    field_ref = angular_spectrum(input_field * np.exp(1j * phase_ref), wavelength, dx, z)
    field_test = angular_spectrum(input_field * np.exp(1j * phase_test), wavelength, dx, z)
    intensity_ref = np.abs(field_ref) ** 2
    intensity_test = np.abs(field_test) ** 2
        # Після поширення обчислюємо показники, якими кількісно порівнюємо два поля.
    metrics = calculate_index(field_ref, field_test, roi_fraction)

    return {"x_um": coords * 1e6, "phase_ref": phase_ref, "phase_test": phase_test,
            "field_ref": field_ref, "field_test": field_test,
            "intensity_ref": intensity_ref, "intensity_test": intensity_test,
                        # Для відображення кожна карта нормується до власного максимуму; це лише візуалізація.
            "display_intensity_ref": intensity_ref / intensity_ref.max() if intensity_ref.max() > 0 else intensity_ref,
            "display_intensity_test": intensity_test / intensity_test.max() if intensity_test.max() > 0 else intensity_test,
            "metrics": metrics, "kind": kind, "levels": levels,
            "wavelength_nm": wavelength_nm, "pixel_um": pixel_um,
            "size": size, "focal_mm": focal_mm, "z_mm": z_mm,
            "roi_fraction": roi_fraction}


def run_quantization_sweep(output_csv: str | Path, **kwargs) -> list[dict[str, float]]:
    """Обчислює показники для L=2,4,...,256 і зберігає CSV."""
    levels_list = [2, 4, 8, 16, 32, 64, 128, 256]
    rows = []
    for levels in levels_list:
        result = simulate_pair(levels=levels, **kwargs)
        row = {"levels": levels, **result["metrics"]}
        rows.append(row)
    path = Path(output_csv)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return rows


# Цей клас створює графічне вікно, поля введення, кнопки, зображення та блок результатів.
class MetaSurfaceApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title(APP_NAME)
        self.root.geometry("1360x900")
        self.root.minsize(1050, 720)
        self.result = None
        self.photos = {}
        self._build()

    def _entry(self, parent, label, default, width=8):
        ttk.Label(parent, text=label).pack(side="left", padx=(7, 2))
        var = tk.StringVar(value=str(default))
        ttk.Entry(parent, textvariable=var, width=width).pack(side="left", padx=(0, 4))
        return var

    def _build(self):
        controls = ttk.LabelFrame(self.root, text="Параметри моделі", padding=8)
        controls.pack(fill="x", padx=8, pady=6)
        row1 = ttk.Frame(controls); row1.pack(fill="x", pady=3)
        ttk.Label(row1, text="Профіль:").pack(side="left")
        self.kind = tk.StringVar(value="Лінза")
        ttk.Combobox(row1, textvariable=self.kind, state="readonly", width=20,
                     values=["Лінза", "Аксикон", "Вихрова пластинка", "Керований нахил", "Випадкова фаза"]).pack(side="left", padx=4)
        self.wavelength = self._entry(row1, "λ, нм:", 532)
        self.pixel = self._entry(row1, "Крок, мкм:", 4)
        self.size = self._entry(row1, "Сітка N×N:", 256)
        self.focal = self._entry(row1, "Фокус, мм:", 20)
        self.distance = self._entry(row1, "Відстань z, мм:", 20)
        row2 = ttk.Frame(controls); row2.pack(fill="x", pady=3)
        self.beam = self._entry(row2, "Радіус пучка, мкм:", 150)
        self.levels = self._entry(row2, "Фазові рівні L:", 8)
        self.angle = self._entry(row2, "Кут аксикона, °:", 1.0)
        self.charge = self._entry(row2, "Заряд вихору m:", 1)
        self.random = self._entry(row2, "Множник фази:", 1.0)
        self.seed = self._entry(row2, "Seed:", 2026)
        self.aperture = tk.BooleanVar(value=True)
        ttk.Checkbutton(row2, text="Кінцева апертура", variable=self.aperture).pack(side="left", padx=6)
        row3 = ttk.Frame(controls); row3.pack(fill="x", pady=4)
        ttk.Button(row3, text="Порівняти поля та обчислити ІБВП", command=self.run).pack(side="left", padx=3)
        ttk.Button(row3, text="Експеримент L=2…256 → CSV", command=self.export_sweep).pack(side="left", padx=3)
        ttk.Button(row3, text="Зберегти карти PNG", command=self.save_maps).pack(side="left", padx=3)
        ttk.Label(row3, text="Скалярна модель; не повний розв'язок Максвелла.", foreground="#8a3d00").pack(side="right")

        views = ttk.Frame(self.root, padding=6); views.pack(fill="both", expand=True)
        self.view_labels = {}
        for key, title in [("phase", "Квантизований фазовий профіль"),
                           ("ref", "Еталонна інтенсивність"),
                           ("test", "Інтенсивність після квантування")]:
            frame = ttk.LabelFrame(views, text=title, padding=5)
            frame.pack(side="left", fill="both", expand=True, padx=4)
            label = ttk.Label(frame, text="Запустіть розрахунок", anchor="center")
            label.pack(fill="both", expand=True)
            self.view_labels[key] = label

        out_frame = ttk.LabelFrame(self.root, text="Метрики та результати", padding=7)
        out_frame.pack(fill="x", padx=8, pady=6)
        self.output = tk.Text(out_frame, height=8, font=("Consolas", 10), wrap="word")
        self.output.pack(fill="x")
        self.output.insert("1.0", "Після запуску тут з'являться E_i, E_phi, E_eta та ІБВП.")
        self.output.configure(state="disabled")

            # Зчитуємо введені користувачем значення й перетворюємо їх на числа для розрахунку.
    def _read_params(self):
        return dict(kind=self.kind.get(), wavelength_nm=float(self.wavelength.get()),
                    pixel_um=float(self.pixel.get()), size=int(self.size.get()),
                    focal_mm=float(self.focal.get()), z_mm=float(self.distance.get()),
                    beam_radius_um=float(self.beam.get()), levels=int(self.levels.get()),
                    axicon_angle_deg=float(self.angle.get()), vortex_charge=int(self.charge.get()),
                    random_strength=float(self.random.get()), seed=int(self.seed.get()),
                    aperture_enabled=bool(self.aperture.get()))

    @staticmethod
            # Перетворюємо числовий масив на зображення для показу у вікні програми.
    def _image(array, phase=False):
        arr = np.asarray(array, dtype=float)
        if phase:
            norm = np.clip((arr + np.pi) / (2 * np.pi), 0, 1)
        else:
                        # Логарифмічне відображення допомагає побачити слабкі деталі дифракційної картини.
            # log display improves visibility of weaker diffraction structure
            norm = np.log1p(99 * np.clip(arr, 0, 1)) / np.log(100)
        img = Image.fromarray((norm * 255).astype(np.uint8), mode="L")
        img.thumbnail((400, 390), Image.Resampling.NEAREST)
        return ImageTk.PhotoImage(img.convert("RGB"))

            # Виконуємо один розрахунок і показуємо фазовий профіль, обидві карти та метрики.
    def run(self):
        try:
            self.result = simulate_pair(**self._read_params())
            r = self.result
            imgs = {"phase": self._image(r["phase_test"], True),
                    "ref": self._image(r["display_intensity_ref"]),
                    "test": self._image(r["display_intensity_test"])}
            self.photos = imgs
            for key, photo in imgs.items():
                self.view_labels[key].configure(image=photo, text="")
            m = r["metrics"]
            lines = [
                f"Профіль: {r['kind']} | L={r['levels']} | λ={r['wavelength_nm']:g} нм | dx={r['pixel_um']:g} мкм | сітка {r['size']}×{r['size']}",
                f"E_i   — нормована похибка інтенсивності:             {m['E_i']:.8f}",
                f"E_phi — зважена фазова неузгодженість:              {m['E_phi']:.8f}",
                f"E_eta — різниця часток енергії в центральній ROI:   {m['E_eta']:.8f}",
                f"ІБВП D_B = 1-(1-E_i)(1-E_phi)(1-E_eta):             {m['D_B']:.8f}",
                f"Частка енергії в ROI: еталон={m['eta_ref']:.6f}, тест={m['eta_test']:.6f}",
            ]
            self.output.configure(state="normal"); self.output.delete("1.0", "end")
            self.output.insert("1.0", "\n".join(lines)); self.output.configure(state="disabled")
        except Exception as exc:
            messagebox.showerror("Помилка розрахунку", str(exc))

            # Запускаємо серію для L=2,4,...,256 і записуємо результати у вибраний CSV-файл.
    def export_sweep(self):
        try:
            params = self._read_params()
            params.pop("levels", None)
            path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")],
                                                initialfile="ibvp_quantization_sweep.csv",
                                                title="Зберегти результати експерименту")
            if not path:
                return
            rows = run_quantization_sweep(path, **params)
            messagebox.showinfo("Готово", f"Збережено {len(rows)} сценаріїв у файлі:\n{path}")
        except Exception as exc:
            messagebox.showerror("Помилка експерименту", str(exc))

            # Зберігаємо три зображення на диск, щоб їх можна було вставити у звіт.
    def save_maps(self):
        if self.result is None:
            messagebox.showinfo("Немає результату", "Спочатку запустіть розрахунок.")
            return
        directory = filedialog.askdirectory(title="Виберіть папку для PNG-карт")
        if not directory:
            return
        r = self.result
        try:
            from PIL import Image
            maps = [("phase_quantized.png", (r["phase_test"] + np.pi) / (2*np.pi)),
                    ("intensity_reference.png", r["display_intensity_ref"]),
                    ("intensity_quantized.png", r["display_intensity_test"])]
            for name, arr in maps:
                Image.fromarray((np.clip(arr, 0, 1) * 65535).astype(np.uint16)).save(Path(directory) / name)
            messagebox.showinfo("Готово", f"Карти збережено в:\n{directory}")
        except Exception as exc:
            messagebox.showerror("Помилка збереження", str(exc))


# Автоматичні перевірки контролюють базові властивості формул і повторюваність розрахунку.
def _self_test():
    rng = np.random.default_rng(123)
    u = rng.normal(size=(32, 32)) + 1j * rng.normal(size=(32, 32))
        # Якщо поля однакові, усі відхилення мають бути практично нульовими.
    same = calculate_index(u, u)
    assert max(same["E_i"], same["E_phi"], same["E_eta"], same["D_B"]) < 1e-12, same
        # Спільний фазовий зсув не повинен створювати фазову похибку за нашим визначенням.
    global_phase = calculate_index(u, u * np.exp(1j * 1.234))
    assert global_phase["E_phi"] < 1e-12 and global_phase["D_B"] < 1e-12, global_phase
    zero = np.zeros((8, 8), dtype=np.complex128)
    assert calculate_index(zero, zero)["D_B"] == 0.0
    assert calculate_index(zero, u[:8, :8])["E_phi"] == 1.0
        # Перевіряємо, що однакові параметри й seed дають однаковий результат.
    # Determinism check on a small, fixed scenario.
    a = simulate_pair(size=128, levels=8, seed=2026)
    b = simulate_pair(size=128, levels=8, seed=2026)
    assert abs(a["metrics"]["D_B"] - b["metrics"]["D_B"]) < 1e-14
    print("Усі автоматичні перевірки пройдено.")


# Точка входу: програма або відкриває графічний інтерфейс, або запускає тест/експорт CSV.
def main():
    parser = argparse.ArgumentParser(description="MetaSurface Lab із розрахунком ІБВП")
    parser.add_argument("--self-test", action="store_true", help="запустити математичні перевірки")
    parser.add_argument("--sweep", metavar="CSV", help="виконати експеримент L=2..256 і зберегти CSV")
    parser.add_argument("--size", type=int, default=256, help="розмір сітки для --sweep")
    args = parser.parse_args()
    if args.self_test:
        _self_test(); return
    if args.sweep:
        rows = run_quantization_sweep(args.sweep, size=args.size)
        for row in rows:
            print(row)
        print(f"Результати збережено: {args.sweep}")
        return
    root = tk.Tk()
    try:
        style = ttk.Style()
        if "vista" in style.theme_names(): style.theme_use("vista")
    except tk.TclError:
        pass
    MetaSurfaceApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
