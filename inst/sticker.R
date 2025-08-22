library(hexSticker)
library(showtext)
font_add_google("Lato")

showtext_opts(dpi = 600)

sticker(
  "inst/projection_simple.png",
  s_x = 1,
  s_y = 0.73,
  # s_width=2, s_height=1.5,
  s_width = 0.45,
  # s_height = 0.1,
  package = "cyDefine",
  p_color = "#cd00d0",
  p_family = "Lato",
  p_fontface = "bold",
  p_size = 8,
  h_size = .8,
  spotlight = T,
  l_y = 1.3,
  l_x = 1,
  l_width = 4,
  l_height = 4,
  l_alpha = .3,
  # p_y = 1.5,
  h_fill = "#ffdaba",
  h_color = "#f67b00",
  filename = "inst/cyDefine.png",
  dpi = 600
)

