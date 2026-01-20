# create hex sticker for package
library(hexSticker)
library(here)
hex_sticker_path <- here("package_logo", "hex_logo.png")
imgurl <- here("package_logo", "brech_logo.png")
sticker(imgurl,
        package = "brech",
        p_y = .4, p_size = 24, p_color = "#1A0757",
        s_x = 0.95, s_y = 1.0, s_width = .6, s_height = 0.45,
        h_fill = "#FFFFFF", h_color = "#003E79",
        filename = hex_sticker_path
)
usethis::use_logo(hex_sticker_path)
