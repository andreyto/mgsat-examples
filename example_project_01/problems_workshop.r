
## location of MGSAT code
MGSAT_SRC = "~/work/mgsat"
source(paste(MGSAT_SRC,"dependencies.r",sep="/"),local=T)
## loads dependency packages (which already must be installed)
load_required_packages()

## loads MGSAT code
source(paste(MGSAT_SRC,"report_pandoc.r",sep="/"),local=T)
source(paste(MGSAT_SRC,"power_and_tests.r",sep="/"),local=T)

library(fitdistrplus)

m_a = read.table.m_a("data/1.1.1.1.1-13d859a75539samples.raw.16s.l.2")

#method=c("ident","norm.prop","norm.ihs.prop","norm.clr","norm.boxcox")
m_a = norm.count.m_a(m_a,method="norm.boxcox")

x = m_a$count[,"Bacteroidetes"]
y = m_a$attr[,"DietStatus"]
print(show.distr.group(x = x, group = y))
x.group = x[y=="before.diet"]
descdist(x.group,boot=1000)
plotdist(x.group, "norm", para=list(mean=mean(x.group), sd=sd(x.group)))