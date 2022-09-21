lean <- read.table("data.txt", header = TRUE)
colnames(lean)[17] <- 'race'


temp <- data.frame(table(lean$ID))
temp$Var1 <- as.numeric(as.character(temp$Var1))



# 找出baseline，4个月和9个月都有测量的数据
index <- temp$Var1[temp$Freq == 3]
lean_not_missing <- lean[lean$ID %in% index, ]

# 计算身高和体重三次测量的均值
lean_not_missing$Height <- apply(lean_not_missing[, 8:10], 1, mean)
lean_not_missing$Weight <- apply(lean_not_missing[, 11:13], 1, mean)

# 计算SBP和DBP三次测量的均值
lean_not_missing$SBP_mean <- apply(lean_not_missing[, 2:4], 1, mean)
lean_not_missing$DBP_mean <- apply(lean_not_missing[, 5:7], 1, mean)


# 计算BMI指数
lean_not_missing$BMI <- lean_not_missing$Weight/(lean_not_missing$Height^2) * 703


# 复制lean_not_missing
dat <- lean_not_missing

# 去掉第二次测量的收缩压和舒张压，身高和体重
dat <- dat[, -c(3, 6)]
dat <- dat[, -c(6:11)]
dat <- dat[, -c(11,12)]


# 计算第9个月与baseline的BMI变化率
dat$ID <- factor(dat$ID)

BmiChange <- function(x) {
  (x[3] - x[1])/x[1]
}

BMI_change <- tapply(dat$BMI, dat$ID, BmiChange)


# 取出dat中表示第9个月的数据

nine_index <- seq(3, 369, by = 3)

nine_dat <- dat[nine_index, ]

nine_dat$BMI_change <- BMI_change


nine_dat$ID <- 1:123
# 去掉第80个人
nine_dat <- nine_dat[-80, ]
## 再去掉第33， 58个人
nine_dat <- nine_dat[-c(33, 57), ]


nine_dat$Group[nine_dat$Group %in% c(1)] = 0
nine_dat$Group[nine_dat$Group %in% c(2, 3, 4)] = 1


remove(lean)

save(nine_dat, file = 'nine_dat.RData')

