calc.probs<-function(data,interval=1) {	cat('Are you calculating probabilities for:\n')	cat('(1) Regression models or (2) Variance component models?\n')	choice <- scan(n = 1)		if (choice == 1) {		require(cnf2freq)
	raw<- .External("cnf2doer", data, resolution=interval, PACKAGE = 'cnf2freq')
	n<-dim(raw[[1]])[1]/4
	ch<-sort(unique(data$geno$chr))
	id<-sapply(rownames(raw[[1]][seq(1,4*n-3,by=4),]),function(x) strsplit(x,"\\.")[[1]][1],USE.NAMES=FALSE)
	orig.ord<-rownames(data$pheno)[which(data$pheno$generation==3)]
	reorder<-unlist(sapply(orig.ord,function(x,y) {which(y==x)},y=id,USE.NAMES=F))
	vals<-lapply(raw,calc.chrom.probs,n=n,id=id,s=interval,order=reorder)
	names(vals)<-paste("chrom.",ch,sep="")
	if (!is.na(data$sex.chrom)) { 
		sex.chr<-vals[[which(ch==data$sex.chrom)]]
		if (length(ch)>1) {
			auto<-vals[-which(ch==data$sex.chrom)]
			auto.names<-names(vals)[-which(ch==data$sex.chrom)]
			exist.auto=1
		} else {
			exist.auto<-0
		}
	} else {
		auto<-vals
		exist.auto<-1
	}
	if (data$backcross==1) {
		bc.origin<-sapply(orig.ord,backcross.origin,pheno=data$pheno,heterogam.sex=data$heterogam)
		if (exist.auto==1) {
			errors<-sapply(auto,bc.auto.check,origin=bc.origin,n=n)
			summary.err<-rowSums(errors)
			if (length(which(summary.err>0))>1) {
				err<-errors[which(summary.err>0),]
				print(paste(dim(err)[1],"individuals have non-zero probabilities that do not match with the pedigree"))
				data$prob.errors<-err
			}		
			data$ad<-lapply(auto,create.bc.ad,backcross.line=data$backcross.line)
			data$origin.probs<-auto
		}
		if (!is.na(data$sex.chrom)) {
			data$sex.chrom.origin.probs<-sex.chr
			if (!is.na(data$backcross.parent) & data$backcross.parent!=data$heterogam) {
				print("Warning, no recombination between lines for the sex chromosome so cannot map QTL on sex chromosome")
				if (!is.na(data$backcross.line)) {
					if (data$sex.restrict==1) {
						print ("Warning all offspring of same sex inherit same line origin of sex chromosome(s) so cannot estimate sex chromosome effects")
					} else {
						print("Warning all heterozygous offspring have same line origin of sex chromosome so cannot estimate sex chromosome effects in this sex")
					}
				}
			}
			if (is.na(data$backcross.parent) | data$backcross.parent==data$heterogam | data$sex.restrict!=1 | is.na(data$backcross.line)) {
				p11<-rowSums(sex.chr$prob.11)
				p12<-rowSums(sex.chr$prob.12)
				p21<-rowSums(sex.chr$prob.21)
				p22<-rowSums(sex.chr$prob.22)
				sex.ad<-sapply(1:n,create.bc.sex.ad,probs=sex.chr,p11=p11,p12=p12,p21=p21,p22=p22,origin=bc.origin,heterogam.sex=data$heterogam,sex.restrict=data$sex.restrict,backcross.parent=data$backcross.parent,backcross.line=data$backcross.line)
				if (is.na(data$backcross.parent) | data$backcross.parent==data$heterogam) {
					if (dim(sex.ad)[1]==4) {
						a1<-data.frame(matrix(unlist(sex.ad[1,]),byrow=TRUE,nrow=n))
						a2<-data.frame(matrix(unlist(sex.ad[2,]),byrow=TRUE,nrow=n))
						d2<-data.frame(matrix(unlist(sex.ad[3,]),byrow=TRUE,nrow=n))
						e<-unlist(sex.ad[4,])
						tot.e<-length(which(e==1))
						if (tot.e>0) {
							print(paste(tot.e,"individuals have non-zero probabilities on the sex chromosome that do not match with the pedigree"))
							names(e)<-orig.ord
							data$sex.prob.errors<-e
						}	
						rownames(a1)<-orig.ord
						colnames(a1)<-paste(seq(0,by=interval,length.out=dim(a1)[2]),"cM")
						rownames(a2)<-orig.ord
						colnames(a2)<-paste(seq(0,by=interval,length.out=dim(a1)[2]),"cM")
						rownames(d2)<-orig.ord
						colnames(d2)<-paste(seq(0,by=interval,length.out=dim(a1)[2]),"cM")
						data$sex.chrom.ad<-list(a1=a1,a2=a2,d2=d2)	
					} else {
						a1<-data.frame(matrix(unlist(sex.ad[1,]),byrow=TRUE,nrow=n))
						a2<-data.frame(matrix(unlist(sex.ad[2,]),byrow=TRUE,nrow=n))
						e<-unlist(sex.ad[3,])
						tot.e<-length(which(e==1))
						if (tot.e>0) {
							print(paste(tot.e,"individuals have non-zero probabilities on the sex chromosome that do not match with the pedigree"))
							names(e)<-orig.ord
							data$sex.prob.errors<-e	
						}
						rownames(a1)<-orig.ord
						colnames(a1)<-paste(seq(0,by=interval,length.out=dim(a1)[2]),"cM")
						rownames(a2)<-orig.ord
						colnames(a2)<-paste(seq(0,by=interval,length.out=dim(a1)[2]),"cM")
						data$sex.chrom.ad<-list(a1=a1,a2=a2)		
					}			
				} else {
					if (dim(sex.ad)[1]==4) {
						a1<-unlist(sex.ad[1,])
						a2<-unlist(sex.ad[2,])
						d2<-unlist(sex.ad[3,])
						e<-unlist(sex.ad[4,])
						tot.e<-length(which(e==1))
						if (tot.e>0) {
							print(paste(tot.e,"individuals have non-zero probabilities on the sex chromosome that do not match with the pedigree"))
							names(e)<-orig.ord
							data$sex.prob.errors<-e
						}
						names(a1)<-orig.ord
						names(a2)<-orig.ord
						names(d2)<-orig.ord	
						data$sex.chrom.ad<-list(a1=a1,a2=a2,d2=d2)	
					} else {
						if (dim(sex.ad)[1]==3) {
							a1<-unlist(sex.ad[1,])
							a2<-unlist(sex.ad[2,])
							e<-unlist(sex.ad[3,])
							names(a1)<-orig.ord
							names(a2)<-orig.ord
							tot.e<-length(which(e==1))
							if (tot.e>0) {
								print(paste(tot.e,"individuals have non-zero probabilities on the sex chromosome that do not match with the pedigree"))
								names(e)<-orig.ord
								data$sex.prob.errors<-e
							}			
							data$sex.chrom.ad<-list(a1=a1,a2=a2)
						} else {
							a2<-unlist(sex.ad[1,])
							e<-unlist(sex.ad[2,])
							names(a2)<-orig.ord
							tot.e<-length(which(e==1))
							if (tot.e>0) {
								print(paste(tot.e,"individuals have non-zero probabilities on the sex chromosome that do not match with the pedigree"))
								names(e)<-orig.ord
								data$sex.prob.errors<-e
							}			
							data$sex.chrom.ad<-list(a2=a2)
						}
					}
				}
			}
		}
	} else {
		if (exist.auto==1) {
			data$ad<-lapply(auto,create.ad)
			data$origin.probs<-auto
		}
		if (!is.na(data$sex.chrom)) {
			sex.origin<-sapply(orig.ord,sex.chrom.origin,pheno=data$pheno,heterogam.sex=data$heterogam)
			p11<-rowSums(sex.chr$prob.11)
			p12<-rowSums(sex.chr$prob.12)
			p21<-rowSums(sex.chr$prob.21)
			p22<-rowSums(sex.chr$prob.22)
			sex.ad<-sapply(1:n,create.sex.ad,probs=sex.chr,p11=p11,p12=p12,p21=p21,p22=p22,origin=sex.origin,heterogam.sex=data$heterogam,sex.restrict=data$sex.restrict)
			a1<-data.frame(matrix(unlist(sex.ad[1,]),byrow=TRUE,nrow=n))
			rownames(a1)<-orig.ord
			colnames(a1)<-paste(seq(0,by=interval,length.out=dim(a1)[2]),"cM")
			a2<-data.frame(matrix(unlist(sex.ad[2,]),byrow=TRUE,nrow=n))
			rownames(a2)<-orig.ord
			colnames(a2)<-paste(seq(0,by=interval,length.out=dim(a2)[2]),"cM")
			if (data$sex.restrict==1) {
				e<-unlist(sex.ad[3,])
				tot.e<-length(which(e==1))
				if (tot.e>0) {
					print(paste(tot.e,"individuals have non-zero probabilities on the sex chromosome that do not match with the pedigree"))
					names(e)<-orig.ord
					data$sex.prob.errors<-e	
				}
				data$sex.chrom.ad<-list(a1=a1,a2=a2)
			} else {
				e<-unlist(sex.ad[4,])
				tot.e<-length(which(e==1))
				if (tot.e>0) {
					print(paste(tot.e,"individuals have non-zero probabilities on the sex chromosome that do not match with the pedigree"))
					names(e)<-orig.ord
					data$sex.prob.errors<-e	
				}
				d2<-data.frame(matrix(unlist(sex.ad[3,]),byrow=TRUE,nrow=n))
				rownames(d2)<-orig.ord
				colnames(d2)<-paste(seq(0,by=interval,length.out=dim(d2)[2]),"cM")
				data$sex.chrom.ad<-list(a1=a1,a2=a2,d2=d2)
			}
			data$sex.chrom.origin.probs<-sex.chr
		}
	}
	data		}		else {		require(cnf2freqibdtracer)		output <<- .External("cnf2doer", data, resolution = interval, PACKAGE = 'cnf2freqibdtracer')	}
}


create.sex.ad<-function(i,p11,p12,p21,p22,probs,origin,heterogam.sex,sex.restrict) {
	x<-dim(probs$prob.21)[2]	
	if(sex.restrict==1) {
		if (heterogam.sex==1) {
			if (origin[1,i]==1) {
				if ((p21[i]+p22[i])>0) {
					e<-1
				} else {
					e<-0
				}
			} else {
				if ((p11[i]+p12[i])>0) {
					e<-1
				} else {
					e<-0
				}
			}
			if (origin[2,i]==1) {
				a1<-probs$prob.11[i,]+probs$prob.21[i,]
				a2<-rep(0,x)
			} else {
				a2<-probs$prob.11[i,]+probs$prob.21[i,]
				a1<-rep(0,x)
			}
		} else {
			if (origin[1,i]==1) {
				if ((p12[i]+p22[i])>0) {
					e<-1
				} else {
					e<-0
				}
			} else {
				if ((p11[i]+p21[i])>0) {
					e<-1
				} else {
					e<-0
				}
			}
			if (origin[2,i]==2) {
				a1<-probs$prob.11[i,]+probs$prob.12[i,]
				a2<-rep(0,x)
			} else {
				a2<-probs$prob.11[i,]+probs$prob.12[i,]
				a1<-rep(0,x)
			}
		}
		list(a1=a1,a2=a2,e=e)
	} else {
		if (heterogam.sex==1) {
			if (origin[1,i]==1) {
				if ((p21[i]+p22[i])>0) {
					e<-1
				} else {
					e<-0
				}
			} else {
				if ((p11[i]+p12[i])>0) {
					e<-1
				} else {
					e<-0
				}
			}
			if (origin[2,i]==1) {
				a1<-probs$prob.11[i,]+probs$prob.21[i,]
				a2<-rep(0,x)
				d1<-rep(0,x)
			} else {
				a2<-probs$prob.11[i,]-probs$prob.22[i,]
				d2<-probs$prob.12[i,]+probs$prob.21[i,]
				a1<-rep(x)
			}
		} else {
			if (origin[1,i]==1) {
				if ((p12[i]+p22[i])>0) {
					e<-1
				} else {
					e<-0
				}
			} else {
				if ((p11[i]+p21[i])>0) {
					e<-1
				} else {
					e<-0
				}
			}
			if (origin[2,i]==2) {
				a1<-probs$prob.11[i,]+probs$prob.12[i,]
				a2<-rep(0,x)
				d2<-rep(0,x)
			} else {
				a2<-probs$prob.11[i,]-probs$prob.22[i,]
				d2<-probs$prob.12[i,]+probs$prob.21[i,]
				a1<-rep(0,x)
			}
		}
		list(a1=a1,a2=a2,d2=d2,e=e)
	}
}


create.ad<-function(probs) {
	a<-data.frame(probs$prob.11-probs$prob.22)
	d<-data.frame(probs$prob.12+probs$prob.21)
	colnames(a)<-colnames(probs$prob.11)
	colnames(d)<-colnames(probs$prob.11)	
	list(a=a,d=d)
}


calc.chrom.probs<-function(ch.raw,n,id,s,order) {
	prob.11<-data.frame(ch.raw[seq(1,4*n-3,by=4),])
	rownames(prob.11)<-id
	l<-dim(prob.11)[2]	
	colnames(prob.11)<-paste(seq(0,by=s,length.out=l),"cM")
	prob.11<-prob.11[order,]
	prob.12<-data.frame(ch.raw[seq(2,4*n-2,by=4),]) 
	rownames(prob.12)<-id
	colnames(prob.12)<-paste(seq(0,by=s,length.out=l),"cM")
	prob.12<-prob.12[order,]
	prob.21<-data.frame(ch.raw[seq(3,4*n-1,by=4),])
	rownames(prob.21)<-id
	colnames(prob.21)<-paste(seq(0,by=s,length.out=l),"cM")
	prob.21<-prob.21[order,]
	prob.22<-data.frame(ch.raw[seq(4,4*n,by=4),])
	rownames(prob.22)<-id
	colnames(prob.22)<-paste(seq(0,by=s,length.out=l),"cM")
	prob.22<-prob.22[order,]
	list(prob.11=prob.11,prob.12=prob.12,prob.21=prob.21,prob.22=prob.22)
}


sex.chrom.origin<-function(pheno,id,heterogam.sex) {
	i<-which(rownames(pheno)==id)
	if (heterogam.sex==1) {
		row.parent<-which(rownames(pheno)==pheno$parent_1[i])
		if (pheno$sex[i]==1) {
			row.grandparent<-which(rownames(pheno)==pheno$parent_1[row.parent])
		} else {
			row.grandparent<-which(rownames(pheno)==pheno$parent_2[row.parent])
		}
	} else {
		row.parent<-which(rownames(pheno)==pheno$parent_2[i])
		if (pheno$sex[i]==2) {
			row.grandparent<-which(rownames(pheno)==pheno$parent_2[row.parent])
		} else {
			row.grandparent<-which(rownames(pheno)==pheno$parent_1[row.parent])
		}
	}
	ori<-pheno$line[row.grandparent]
	sex<-pheno$sex[i]
	list(origin=ori,sex=sex)
}


backcross.origin<-function(pheno,id,heterogam.sex) {
	i<-which(rownames(pheno)==id)
	row.parent.1<-which(rownames(pheno)==pheno$parent_1[i])
	row.parent.2<-which(rownames(pheno)==pheno$parent_2[i])
	row.gp.11<-which(rownames(pheno)==pheno$parent_1[row.parent.1])
	row.gp.12<-which(rownames(pheno)==pheno$parent_2[row.parent.1])
	row.gp.21<-which(rownames(pheno)==pheno$parent_1[row.parent.2])
	row.gp.22<-which(rownames(pheno)==pheno$parent_2[row.parent.2])
	if (pheno$line[row.gp.11]==pheno$line[row.gp.12]) {
		bc.parent=1
		bc.line=pheno$line[row.gp.11]
	} else {
		if (pheno$line[row.gp.21]==pheno$line[row.gp.22]) {
			bc.parent=2
			bc.line=pheno$line[row.gp.21]
		} else {
			print(paste(id,"is not from a backcross"))
		}
	}
	if (heterogam.sex==1) {
		if(pheno$sex[i]==1) {
			row.grandparent<-row.gp.11
		} else {
			row.grandparent<-row.gp.12
		}
	} else {
		if(pheno$sex[i]==1) {
			row.grandparent<-row.gp.21
		} else {
			row.grandparent<-row.gp.22
		}
	}
	ori<-pheno$line[row.grandparent]
	sex<-pheno$sex[i]
	list(bc.parent=bc.parent,bc.line=bc.line,het.sex.origin=ori,sex=sex)
}
	
bc.auto.check<-function(probs,origin,n) {
	p11<-rowSums(probs$prob.11)
	p12<-rowSums(probs$prob.12)
	p21<-rowSums(probs$prob.21)
	p22<-rowSums(probs$prob.22)
	sapply(1:n,bc.prob.check,p11=p11,p12=p12,p21=p21,p22=p22,origin=origin)
}

bc.prob.check<-function(i,p11,p12,p21,p22,origin) {
	if (origin[1,i]==1) {
		if (origin[2,i]==1) {
			if ((p21[i]+p22[i])>0) {
				e<-1
			} else {
				e<-0
			}
		} else {
			if ((p11[i]+p12[i])>0) {
				e<-1
			} else {
				e<-0
			}
		}
	} else {
		if (origin[2,i]==1) {
			if ((p12[i]+p22[i])>0) {
				e<-1
			} else {
				e<-0
			}
		} else {
			if ((p11[i]+p21[i])>0) {
				e<-1
			} else {
				e<-0
			}
		}
	}
	e	
}


create.bc.ad<-function(probs,backcross.line) {
	if(is.na(backcross.line)) {
		a<-data.frame(probs$prob.11-probs$prob.22)
		d<-data.frame(probs$prob.12+probs$prob.21)
		colnames(a)<-colnames(probs$prob.11)
		colnames(d)<-colnames(probs$prob.11)
		list(a=a,d=d)
	} else {
		if (backcross.line==1) {
			a<-probs$prob.11
			list(a=a)
		} else	{
			a<-probs$prob.21+probs$prob.12
			list(a=a)
		} 
	}	
}

create.bc.sex.ad<-function(i,p11,p12,p21,p22,probs,origin,heterogam.sex,sex.restrict,backcross.parent,backcross.line) {
	if (origin[1,i]!=heterogam.sex) {	
		if (origin[1,i]==1) {
			if (origin[2,i]==1) {
				if (origin[3,i]==1) {
					if ((p12[i]+p21[i]+p22[i])>0) {
						e<-1
					} else {
						e<-0
					}
				} else {
					if ((p11[i]+p21[i]+p22[i])>0) {
						e<-1
					} else {
						e<-0
					}
				}
			} else {
				if (origin[3,i]==1) {
					if ((p12[i]+p11[i]+p22[i])>0) {
						e<-1
					} else {
						e<-0
					}
				} else {
					if ((p11[i]+p21[i]+p12[i])>0) {
						e<-1
					} else {
						e<-0
					}
				}
			}
		} else {
			if (origin[2,i]==1) {
				if (origin[3,i]==1) {
					if ((p12[i]+p21[i]+p22[i])>0) {
						e<-1
					} else {
						e<-0
					}
				} else {
					if ((p11[i]+p12[i]+p22[i])>0) {
						e<-1
					} else {
						e<-0
					}
				}
			} else {
				if (origin[3,i]==1) {
					if ((p11[i]+p21[i]+p22[i])>0) {
						e<-1
					} else {
						e<-0
					}
				} else {
					if ((p11[i]+p12[i]+p21[i])>0) {
						e<-1
					} else {
						e<-0
					}
				}
			}
		}
	} else {	
		if (origin[1,i]==1) {
			if (origin[2,i]==1) {
				if ((p21[i]+p22[i])>0) {	
					e<-1
				} else {
					e<-0
				}
			} else {
				if ((p11[i]+p12[i])>0) {
					e<-1
				} else {
					e<-0
				}
			}
		} else {
			if (origin[2,i]==1) {
				if ((p12[i]+p22[i])>0) {
					e<-1
				} else {
					e<-0
				}
			} else {
				if ((p11[i]+p21[i])>0) {
					e<-1
				} else {
					e<-0
				}
			}
		}
	}
	if (!is.na(backcross.parent) & backcross.parent!=heterogam.sex) {
		if (sex.restrict==1) {
			if (is.na(backcross.line)) {
				if (origin[4,i]!=heterogam.sex) {
					a1<-0
					if (origin[2,i]<-1) {
						a2<-1
					} else {
						a2<-0
					}
				} else {
					a2<-0
					if (origin[2,i]<-1) {
						a1<-1
					} else {
						a1<-0
					}
				}
				list(a1=a1,a2=a2,e=e)
			}
		} else {
			if (!is.na(backcross.line)) {
				if (origin[4,i]==heterogam.sex) {
					a2<-0
				} else {
					if (origin[3,i]<-1) {
						a2<-1
					} else {
						a2<-0
					}
				}
				list(a2=a2,e=e)
			} else {
				if (origin[4,i]==heterogam.sex) {
					if (origin[2,i]<-1) {
						a1<-1
						a2<-0
						d2<-0
					} else {
						a1<-0
						a2<-0
						d2<-0
					}
				} else {
					if (origin[2,i]<-1) {
						if (origin[3,i]<-1) {
							a2<-1
							d2<-0
						} else {
							d2<-1
							a2<-0
						}
					} else {
						if (origin[3,i]<-1) {
							d2<-1
							a2<-0
						} else {
							a2<--1
							d2<-0
						}
					}
					a1<-0
				}			
				list(a1=a1,a2=a2,d2=d2,e=e)
			}
		}
	} else {
		x<-dim(probs$prob.21)[2]
		if (!is.na(backcross.line)) {
			if (origin[4,i]==heterogam.sex) {
				if (heterogam.sex==1) {	
					a1<-probs$prob.11[i,]+probs$prob.21[i,]		
					a2<-rep(0,x)
				} else {
					a1<-probs$prob.11[i,]+probs$prob.12[i,]
					a2<-rep(0,x)
				}
			} else {
				if (origin[1,i]==1) {
					a2<-probs$prob.11[i,]+probs$prob.21[i,]
					a1<-rep(0,x)
				} else {
					a2<-probs$prob.11[i,]+probs$prob.12[i,]
					a1<-rep(0,x)
				}
			}
			list(a1=a1,a2=a2,e=e)
		} else {
			if (origin[4,i]==heterogam.sex) {
				if (heterogam.sex==1) {	
					a1<-probs$prob.11[i,]+probs$prob.21[i,]		
					a2<-rep(0,x)
					d2<-rep(0,x)
				} else {
					a1<-probs$prob.11[i,]+probs$prob.12[i,]
					a2<-rep(0,x)
					d2<-rep(0,x)
				}
			} else {		
				a2<-probs$prob.11[i,]-probs$prob.22[i,]
				d2<-probs$prob.12[i,]+probs$prob.21[i,]
				a1<-rep(0,x)
			}
			list(a1=a1,a2=a2,d2=d2,e=e)
		}
	}
}

autosome.genome.scan<-function(pheno,data,factors=NULL,covariates=NULL,bk.QTL=NULL,error=normal) { 			data.name<-substitute(data)			data$pheno<-data$pheno[which(data$pheno$generation==3),]			y<-paste("data$pheno$",pheno,sep="")			terms<-""			tms<-c(factors,covariates)			rows<-c("Intercept",factors,covariates)			if (!is.null(factors)) {						for (i in 1:length(factors)) {									x<-which(colnames(data$pheno)==factors[i])									data$pheno[,x]<-as.factor(data$pheno[,x])								}						terms<-paste("data$pheno$",factors,sep="",collapse="+")					}			if (!is.null(covariates)) {						add<-paste("data$pheno$",covariates,sep="",collapse="+")						if (terms=="") {									terms<-add								} else {									terms<-paste(terms,add,sep="+")								}					} 			if (!is.null(bk.QTL)) {						for (i in 1:(length(bk.QTL)/2)) {									inter1<-colnames(data$ad[[1]]$a)									inter2<-strsplit(inter1[2]," ")									interval<-as.numeric(inter2[[1]][1])									pos<-bk.QTL[2*i]/interval+1									tms<-c(tms,paste("QTL",bk.QTL[2*i-1],".",bk.QTL[2*i],sep=""))									if (data$backcross==1 & !is.na(data$backcross.line)) {												rows<-c(rows,paste("QTL",bk.QTL[2*i-1],".",bk.QTL[2*i],".a",sep=""))												if (terms=="") {															terms<-paste("data$ad$chrom.",bk.QTL[2*i-1],"$a[,",pos,"]",sep="")														} else {															terms<-paste(terms,"+data$ad$chrom.",bk.QTL[2*i-1],"$a[,",pos,"]",sep="")														}											} else {												rows<-c(rows,paste("QTL",bk.QTL[2*i-1],".",bk.QTL[2*i],".a",sep=""),paste("QTL",bk.QTL[2*i-1],".",bk.QTL[2*i],".d",sep=""))												if (terms=="") {															terms<-paste("data$ad$chrom.",bk.QTL[2*i-1],"$a[,",pos,"]+data$ad$chrom.",bk.QTL[2*i-1],"$d[,",pos,"]",sep="")														} else {															terms<-paste(terms,"+data$ad$chrom.",bk.QTL[2*i-1],"$a[,",pos,"]+data$ad$chrom.",bk.QTL[2*i-1],"$d[,",pos,"]",sep="")														}											}								}					}			if (terms=="") {						terms<-"1"					} else {						if (length(tms)>1) {									last<-tms[length(tms)]									tms<-tms[-length(tms)]									tms<-paste(paste(tms,collapse=", "),"and",last)								}					}			form<-paste(y,"~",terms,sep="")				red.mod<-lm(formula(form))			if (length(tms)>=1) {						cat ("The results for the starting model fitting",tms,"to",pheno,"in",data.name,"dataset are:\n\n")						out<-signif(summary(red.mod)$coefficients,4)						rownames(out)<-rows						print (out) 						cat("\nF value",signif(summary(red.mod)$fstatistic[1],4),"on",summary(red.mod)$fstatistic[2],"and",summary(red.mod)$fstatistic[3],"degrees of freedom\n\n")					}			red.sse<-sum(red.mod$residuals^2)			red.df<-red.mod$df.residual			res<-lapply(data$ad,norm.chrom.scan,data=data,form=form,red.sse=red.sse,red.df=red.df,rows=rows)			all<-res[[1]]$F.values			ext.pl<-names(res[[1]]$F.values)[length(res[[1]]$F.values)]			sp.pl<-strsplit(ext.pl," ")			pl<-as.numeric(sp.pl[[1]][1])			if (length(res)>1) {						borders<-pl						for (i in 2:length(res)) {									all<-c(all,res[[i]]$F.values)									ext.pl<-names(res[[i]]$F.values)[length(res[[i]]$F.values)]									sp.pl<-strsplit(ext.pl," ")									pl<-pl+1+as.numeric(sp.pl[[1]][1])									borders<-c(borders,pl)								}						max=max(all,na.rm=T)						hi<-numeric()						for (i in 1:length(res)) {									if (res[[i]]$max.F==max) {												hi<-c(hi,i,res[[i]]$pos.max)											}								}						n.hi<-length(hi)/2						if (n.hi==1) {									pos<-hi[2]									mx<-noquote(paste(names(res)[as.numeric(hi[1])],pos))									est<-res[[as.numeric(hi[1])]]$max.QTL.effects								} else {									mx<-""									est<-list()									for (i in 1:(n.hi-1)) {												pos<-hi[2*i]												est[[i]]<-res[[as.numeric(hi[2*i-1])]]$max.QTL.effects												j.p<-strsplit(pos," ")[[1]][1]												names(est)[i]<-paste(names(res)[as.numeric(hi[2*i-1])],j.p,"cM",sep=".")												if (mx=="") {															mx<-paste(names(res)[as.numeric(hi[2*i-1])],pos)														} else {															mx<-paste(mx,", ",names(res)[as.numeric(hi[2*i-1])]," ",pos,sep="")														}											}									pos<-hi[2*n.hi]									mx<-noquote(paste(mx," and ",names(res)[as.numeric(hi[2*n.hi-1])]," ",pos,sep=""))									est[[n.hi]]<-res[[as.numeric(hi[2*n.hi-1])]]$max.QTL.effects									j.p<-strsplit(pos," ")[[1]][1]									names(est)[n.hi]<-paste(names(res)[as.numeric(hi[2*n.hi-1])],j.p,"cM",sep=".")								}						names(all)<-paste(seq(0,pl,length.out=length(all)),"cM")												plot(seq(0,pl,length.out=length(all)),all,typ="l",bty="n",ylab="F value",xlab="position (cM)",main=paste("QTL results for",pheno,"in",data.name,"dataset"))						high<-par("yaxp")[2]						for (i in 1:(length(borders)-1)) {									lines(c(borders[i],borders[i]),c(0,high),lty="dashed")								}						list(max.F.value=max,position.max.F=mx,estimated.effects.max.F=est,variables.in.starting.model=noquote(ifelse(is.null(tms),"none",tms)),all.F.values=all,chrom.boundaries=borders,results.by.chromosome=res)					} else {												plot(seq(0,pl,length.out=length(all)),all,typ="l",bty="n",ylab="F value",xlab="position (cM)",main=paste("QTL results for",pheno,"in",data.name,"dataset"))						list(max.F.value=res[[1]]$max.F,position.max.F=res[[1]]$pos.max,estimated.effects.max.F=res[[1]]$max.QTL.effects,variables.in.starting.model=noquote(ifelse(is.null(tms),"none",tms)),F.values=res[[1]]$F.values)					}		}
norm.chrom.scan<-function(chrom.ad,form,data,red.sse,red.df,rows) {			if (length(chrom.ad)==2) {						res<-mapply(norm.ad.one.pos,chrom.ad$a,chrom.ad$d,MoreArgs=list(form=form,data=data,red.sse=red.sse,red.df=red.df,rows=rows))							rows<-c(rows,"a","d")					} else {						res1<-apply(chrom.ad$a,2,norm.a.one.pos,form=form,data=data,red.sse=red.sse,red.df=red.df,rows=rows)						res<-matrix(unlist(res1,recursive=FALSE),nrow=2)						rows<-c(rows,"a")						colnames(res)<-names(chrom.ad$a)					}			F<-unlist(res[1,])			max<-max(F,na.rm=TRUE)			ind<-which(F==max)			pos<-noquote(names(ind))			est<-res[2,][[ind]]			rownames(est)<-rows			list(F.values=F,max.F=max,pos.max=pos,max.QTL.effects=est)		}norm.ad.one.pos<-function(a,d,form,data,red.sse,red.df,rows) {			new.form<-paste(form,"+a+d",sep="")			model<-lm(formula(new.form))			sse<-sum(model$residuals^2)			df<-model$df.residual			F<-((red.sse-sse)/(red.df-df))/(sse/df)			if(is.nan(F)) {						F<-NA					}			est<-signif(summary(model)$coefficients,4)[,c(1,2,4)]			list(F.values=F,est.effects=est)		}norm.a.one.pos<-function(a,form,data,red.sse,red.df,rows) {			new.form<-paste(form,"+a",sep="")			model<-lm(formula(new.form))			sse<-sum(model$residuals^2)			df<-model$df.residual			F<-((red.sse-sse)/(red.df-df))/(sse/df)			if(is.nan(F)) {						F<-NA					}			est<-signif(summary(model)$coefficients,4)[,c(1,2,4)]			list(F.values=F,est.effects=est)		}
	






