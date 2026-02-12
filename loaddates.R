#load in the date information:

# Lineage 4 data: 210 HIV-, 289 HIV+, 14 unknown
# Lineage 1 data: 27 HIV-, 70 HIV+, 2 unknown
# Lineage 2 data: 10 HIV-, 18 HIV+
# Lineage 3 data: 19 HIV-, 41 HIV+, 1 unknown

load('tms_lin1.Rdata')
load('tms_lin2.Rdata')
load('tms_lin3.Rdata')
load('tms_lin4.Rdata')

colnames(tms.lin1) <- c('date','hiv')
colnames(tms.lin2) <- c('date','hiv')
colnames(tms.lin3) <- c('date','hiv')
colnames(tms.lin4) <- c('date','hiv')

tms.lin1$date <- julian(as.Date(tms.lin1$date))
tms.lin2$date <- julian(as.Date(tms.lin2$date))
tms.lin3$date <- julian(as.Date(tms.lin3$date))
tms.lin4$date <- julian(as.Date(tms.lin4$date))

tms.lin1[tms.lin1$hiv=='Pos','hiv'] <- 'I2'
tms.lin1[tms.lin1$hiv=='Neg','hiv'] <- 'I1'
tms.lin1[tms.lin1$hiv=='Unknown','hiv'] <- 'I1'

tms.lin2[tms.lin2$hiv=='Pos','hiv'] <- 'I2'
tms.lin2[tms.lin2$hiv=='Neg','hiv'] <- 'I1'
tms.lin2[tms.lin2$hiv=='Unknown','hiv'] <- 'I1'

tms.lin3[tms.lin3$hiv=='Pos','hiv'] <- 'I2'
tms.lin3[tms.lin3$hiv=='Neg','hiv'] <- 'I1'
tms.lin3[tms.lin3$hiv=='Unknown','hiv'] <- 'I1'

tms.lin4[tms.lin4$hiv=='Pos','hiv'] <- 'I2'
tms.lin4[tms.lin4$hiv=='Neg','hiv'] <- 'I1'
tms.lin4[tms.lin4$hiv=='Unknown','hiv'] <- 'I1'

