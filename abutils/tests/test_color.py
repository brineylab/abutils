#!usr/env/python
# filename: test_color.py

from ..utils import color


def test_rgb_to_hex():
    black = color.rgb_to_hex((0, 0, 0)).upper()
    red = color.rgb_to_hex((1., 0, 0)).upper()
    green = color.rgb_to_hex((0, 1., 0)).upper()
    blue = color.rgb_to_hex((0, 0, 1.)).upper()
    assert black == '#000000' and red == '#FF0000' and green == '#00FF00' and blue == '#0000FF'


def test_hex_to_rgb():
    black = color.hex_to_rgb('#000000')
    red = color.hex_to_rgb('#FF0000')
    green = color.hex_to_rgb('#00FF00')
    blue = color.hex_to_rgb('#0000FF')
    assert black == (0, 0, 0) and red == (255, 0, 0) and green == (0, 255, 0) and blue == (0, 0, 255)
