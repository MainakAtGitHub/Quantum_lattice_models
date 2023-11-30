function c=clean_BdG_file(filename)
mf=matfile(filename)
mf.Properties.Writable=true;
mf.deltaDiffAcc=0;
mf.deltaMaxAcc=0;
mf.deltaMinAcc=0;
mf.muAcc=0;
acc=mf.nAcc;
acc(length(acc))
mf.nAcc=acc(length(acc));