# 加载必要的包
library(ggplot2)
library(tidyverse)

# 创建数据框
data <- data.frame(
  Orthogroup = c("OG0000000", "OG0000001", "OG0000002", "OG0000003"),
  CPG01 = c(1, 0, 1, 0),
  CPG02 = c(1, 0, 1, 1),
  CPG03 = c(0, 0, 1, 1),
  CPG05 = c(0, 1, 0, 0)
)

# 将数据转换为长格式
data_long <- data %>%
  pivot_longer(cols = -Orthogroup, names_to = "CPG", values_to = "value")

# 按照颜色进行排序
data_long <- data_long %>%
  arrange(value)

# 绘制热图
p <- ggplot(data_long, aes(x = CPG, y = Orthogroup, fill = factor(value))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "blue", "1" = "red")) +
  labs(fill = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Heatmap of Orthogroup vs CPG")

p
ggsave("/bioinfor/mulinlab/wuzhikun/Project/Evolution/orthology/orthoResult/Results_Aug28/Orthogroups/ortho_heatmap.pdf", width=5, height=4)
