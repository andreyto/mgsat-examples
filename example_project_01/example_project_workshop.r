
## Custom metadata loading function (define and pass to proc.project() when default
## implementation load.meta.default() is not sufficient)


load.meta.diet <- function(file.name) {
  ## We are calling default simple table loading here,
  ## but you can do something else. Result must have
  ## rownames() set to SampleID field
    meta = load.meta.default(file.name)
    ## we need to create this as ordered quantile to show
    ## in the right order on plots (and possible used
    ## in models that care about ordered factors)
    meta$age.quant = quantcut.ordered(meta$age)
    meta.aggr = join(meta,
                     ddply(meta,"SubjectID",summarise,
                           visit.4=any(visit==4)),
                     by="SubjectID",
                     match="first")
    stopifnot(!any(is.na(meta.aggr$visit.4)) && 
                nrow(meta.aggr)==nrow(meta))
    
    meta = meta.aggr
    
    return (meta)
}
  

## This function should carry out analysis specific to metadata fields by themselves, without
## relation to the abundance profiles. You can write it to do nothing (empty body).

summary.meta.diet <- function(m_a) {
  
  report$add.header("Summary of metadata variables")
  
  
  report$add.header("Summary of metadata variables after filtering samples")
  
  meta = m_a$attr

  report$add(summary(meta),caption="Summary of metadata variables")
  
  xtabs.formulas = list("~Sample.type+DietStatus","~Drug.Before.Diet + Sample.type",
                        "~Complaints + Sample.type",
                        "~Sample.type+visit","~MatchedGroupID","~Sample.type.1","~SubjectID")
  for(xtabs.formula in xtabs.formulas) {
    fact.xtabs = xtabs(as.formula(xtabs.formula),data=meta,drop.unused.levels=T)
    report$add(fact.xtabs,caption=paste("Sample cross tabulation",xtabs.formula))
    report$add.printed(summary(fact.xtabs))
  }
  
  with(meta,{
    report$add.printed(summary(aov(age~Sample.type)),
                       caption="ANOVA for age and sample type")
    report$add(qplot(Sample.type,age,geom="violin"),
               caption="Violin plot for age and sample type")
  })
  
  with(meta,{
    report$add(cor.test(age,
                                visit,
                                method="spearman"),
                       caption="Spearman RHO for age and visit")
    
  })
  with(meta[meta$Sample.type=="patient",],{
    report$add(cor.test(age,
                                visit,
                                method="spearman"),
                       caption="Spearman RHO for age and visit, patients only")
    
  })
  
  #summary(glht(lm(age~visit,data=meta[meta$Sample.type=="patient",]),linfct="visit=0"))
  #summary(glht(lmer(age~visit+(visit|Sample.type),data=meta),linfct="visit=0"))
  report$add(ggplot(meta,aes(x=visit,y=age,color=Sample.type))+
               geom_point()+
               stat_smooth(method="loess", se = T,degree=1,size=1),
             caption="Plot for age and visit with Loess trend line")
  
}

## This function must generate a list with analysis tasks

gen.tasks.diet <- function() {
  
  task0 = within( mgsat.16s.task.template, {
    
    taxa.levels = c(2)

    descr = "All samples, no aggregation, no tests here, only plots"
    
    read.data.task = within(read.data.task, {
      taxa.summary.file = "example_taxa.summary.file"
      otu.shared.file="example_otu.shared.file"
      cons.taxonomy.file="example_cons.taxonomy.file"
      meta.file="example_meta.txt"
      load.meta.method=load.meta.diet
      load.meta.options=list()    
    })
    
    get.taxa.meta.aggr.base<-function(m_a) { 
      ##any aggregated attributes have to be computed here,
      ##after the available count samples have been joined,
      ##as opposed to in the load.meta() function.
      ## We already have all fields, so not aggregating anything here.
      
      ##Throw away all visits with very few samples
      m_a = subset.m_a(m_a,subset=(m_a$attr$visit<=4))
      return(m_a)
    }
    
    summary.meta.method=summary.meta.diet
    
    test.counts.task = within(test.counts.task, {
      
      norm.count.task = within(norm.count.task, {
        method="norm.ihs.prop"
      })
      
    })
    
  })

  task1 = within( task0, {
    
    main.meta.var = "DietStatus"
    
    descr = "Patients' samples at visits 1 (before diet) and 4 (after diet), only paired samples"
    
    do.summary.meta = F
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$Sample.type=="patient" 
                                   & m_a$attr$visit %in% c(1,4) 
                                   & m_a$attr$visit.1
                                   & m_a$attr$visit.4))
      return(m_a)
    }
    
    test.counts.task = within(test.counts.task, {
      
      #do.divrich = c()
      do.deseq2 = T
      do.adonis = T
      do.genesel = T
      do.stabsel = F
      do.glmer = F
      do.plot.profiles.abund=T
      do.heatmap.abund=T
      do.extra.method = taxa.levels
      
      divrich.task = within(divrich.task,{
        group.attr = main.meta.var
        counts.glm.task = within(counts.glm.task,{
          formula.rhs = main.meta.var
        })
        do.plot.profiles = T
      })
      
      deseq2.task = within(deseq2.task, {
        formula.rhs = main.meta.var
      })
      
      genesel.task = within(genesel.task, {
        group.attr = main.meta.var
        do.plot.profiles = T
        genesel.param = within(genesel.param, {
          block.attr = "SubjectID"
          type="paired"
          replicates=0
        })
      })
      
      
      stabsel.task = within(stabsel.task, {
        resp.attr=main.meta.var
      })
      
      adonis.task = within(adonis.task, {
        
        tasks = list(
          list(formula.rhs=main.meta.var,
               descr="Association with diet status unpaired by subject"),
          list(formula.rhs=main.meta.var,
               strata="SubjectID",
               descr="Association with diet status paired by subject"),
          list(formula.rhs=paste("Drug.Before.Diet * ", main.meta.var),
               strata=NULL,
               descr="Association with Drug use before diet and diet status")
        )
        
        #dist.metr="euclidian"
        #col.trans="standardize"
        
        #norm.count.task = within(norm.count.task, {
        #  method="norm.clr"
        #  drop.features = list()
        #})
        
      })
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c(main.meta.var))
        do.profile=T
        do.clade.meta=F
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c(main.meta.var,"Drug.Before.Diet")
      })
      
      extra.method.task = within(extra.method.task, {
        
        func = function(m_a,m_a.norm,res.tests,norm.count.task.extra) {
          test.dist.matr.within.between(m_a=m_a,
                                        group.attr="DietStatus",
                                        block.attr="SubjectID",
                                        n.perm=8000,
                                        #dist.metr="euclidian",
                                        col.trans="ident",
                                        norm.count.task=norm.count.task.extra
          )
        }
        norm.count.task.extra = within(norm.count.task, {
          method="norm.prop"
          #drop.features = list()
        })
        
      })
      
    })
    
  })
  
  return (list(task1))
}



## number of cores to use on multicore machines
options(mc.cores=4)
options(boot.ncpus=4)
## parallel backend
options(boot.parallel="snow")

## location of MGSAT code
MGSAT_SRC = "~/work/mgsat"

source(paste(MGSAT_SRC,"dependencies.r",sep="/"),local=T)

## Uncomment next line to install packages needed by MGSAT (!!!comment it out
## in all subsequent runs once the packages have been installed!!!).
## Note: you should also pre-install Pandoc program from http://johnmacfarlane.net/pandoc/
## or using your OS package manager (if running on Linux)

#install_required_packages()

## loads dependency packages (which already must be installed)
load_required_packages()

#library("BiocParallel")
#register(SnowParam(4))

## loads MGSAT code
source(paste(MGSAT_SRC,"report_pandoc.r",sep="/"),local=T)
source(paste(MGSAT_SRC,"power_and_tests.r",sep="/"),local=T)

## leave with try.debug=F for production runs
set_trace_options(try.debug=T)

## set incremental.save=T only for debugging or demonstration runs - it forces 
## report generation after adding every header section, thus slowing down
## a long run. But then incremental.save=T, you can open HTML report file in
## a Web browser and refresh it periodically to see it grow.
report <- PandocAT$new(author="noone@mail.com",
                       title="Analysis of Dieting study 16S data",
                       incremental.save=F)


res = proc.project(
  task.generator.method=gen.tasks.diet
)

report$save()

