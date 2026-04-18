"""
Generate prisma_logo.ico (multi-resolution) from prisma_logo.svg.
Run once from the assets/ directory:
    python generate_ico.py
Requires: Pillow  (pip install Pillow)
Optional: cairosvg  (pip install cairosvg)  — for crisp SVG rasterisation.
          Falls back to a pure-Python Pillow drawing if cairosvg is absent.
"""

from pathlib import Path
from PIL import Image
import io
import struct

HERE = Path(__file__).parent
SVG_PATH  = HERE / "prisma_logo.svg"
ICO_PATH  = HERE / "prisma_logo.ico"
PNG_PATH  = HERE / "prisma_logo_256.png"

SIZES = [16, 24, 32, 48, 64, 128, 256]


def render_svg_cairosvg(size: int) -> Image.Image:
    import cairosvg
    png_bytes = cairosvg.svg2png(
        url=str(SVG_PATH),
        output_width=size,
        output_height=size,
    )
    return Image.open(io.BytesIO(png_bytes)).convert("RGBA")


def render_fallback(size: int) -> Image.Image:
    """
    Pure-Pillow fallback: draws the PRISMA prism icon programmatically
    so we never depend on cairosvg being installed.
    Colors match PRISMA v2.0 tokens.
    """
    from PIL import ImageDraw

    img = Image.new("RGBA", (size, size), (4, 7, 13, 255))   # --void
    draw = ImageDraw.Draw(img)
    s = size

    # ── Prism triangle ────────────────────────────────────────────────
    cx = s * 0.352     # apex x
    cy_top = s * 0.109  # apex y
    cy_bot = s * 0.781  # base y
    x_left  = s * 0.109
    x_right = s * 0.594
    triangle = [(cx, cy_top), (x_right, cy_bot), (x_left, cy_bot)]
    draw.polygon(triangle, outline=(255, 255, 255, 50))

    # ── Node ─────────────────────────────────────────────────────────
    nx = cx
    ny = s * 0.445
    nr = max(2, s * 0.016)
    draw.ellipse([nx-nr, ny-nr, nx+nr, ny+nr], fill=(255, 255, 255, 90))

    # ── Beams & dots ─────────────────────────────────────────────────
    channels = [
        (s*0.664, s*0.219, "#7B52FF", 7),
        (s*0.664, s*0.324, "#5BAAFF", 5),
        (s*0.664, s*0.445, "#39FF8A", 9),   # brightest — FITC accent
        (s*0.664, s*0.566, "#FF9B3D", 6),
        (s*0.664, s*0.680, "#FF3D6E", 7),
    ]

    def hex2rgb(h):
        h = h.lstrip("#")
        return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))

    for tx, ty, col, dot_r in channels:
        r, g, b = hex2rgb(col)
        # beam
        draw.line([(nx, ny), (tx, ty)], fill=(r, g, b, 200), width=max(1, s//64))
        # dot
        dr = max(2, s * dot_r / 256)
        draw.ellipse([tx-dr, ty-dr, tx+dr, ty+dr], fill=(r, g, b, 230))

    return img


def get_image(size: int) -> Image.Image:
    try:
        return render_svg_cairosvg(size)
    except Exception:
        return render_fallback(size)


def build_ico(images: dict[int, Image.Image]) -> None:
    """
    Hand-builds a valid ICO file containing PNG-compressed frames
    for every requested size. Works with all sizes including 256.
    """
    num = len(images)
    # ICO header: 6 bytes
    header = struct.pack("<HHH", 0, 1, num)

    # Encode each frame as PNG
    frames = []
    for sz, img in sorted(images.items()):
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        frames.append((sz, buf.getvalue()))

    # Directory entries: 16 bytes each
    dir_offset = 6 + num * 16
    directory = b""
    offset = dir_offset
    for sz, data in frames:
        w = h = sz if sz < 256 else 0   # 0 means 256 in ICO spec
        directory += struct.pack(
            "<BBBBHHII",
            w, h,           # width, height (0 = 256)
            0,              # color count (0 = no palette)
            0,              # reserved
            1,              # planes
            32,             # bit count
            len(data),      # size of image data
            offset,         # offset of image data
        )
        offset += len(data)

    with open(ICO_PATH, "wb") as f:
        f.write(header)
        f.write(directory)
        for _, data in frames:
            f.write(data)


if __name__ == "__main__":
    print("Rendering PRISMA logo at multiple resolutions…")
    images = {}
    for sz in SIZES:
        img = get_image(sz)
        images[sz] = img
        print(f"  {sz}x{sz} OK")

    # Save 256px PNG standalone (for taskbar / about dialog)
    images[256].save(PNG_PATH)
    print(f"Saved {PNG_PATH.name}")

    build_ico(images)
    print(f"Saved {ICO_PATH.name}  ({ICO_PATH.stat().st_size // 1024} KB)")
