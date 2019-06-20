function mtBayesPR_xiolei(genoTrain::DataFrame,genoTrain2::DataFrame, phenoTrain::DataFrame, snpInfo::String, chrs::Int64, fixedRegSize::Int64, varGenotypic::Array{Float64}, varResidual1::Float64,varResidual2::Float64,chainLength::Int64, burnIn::Int64, outputFreq::Int64, onScreen::Bool,resCor::Bool)
    SNPgroups, genoX = prepRegionData(snpInfo, chrs, genoTrain, fixedRegSize)
    these2Keep = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations
    nRegions    = length(SNPgroups)
    println("number of regions: ", nRegions)
    dfEffect    = 4.0
    dfRes       = 4.0
    X1           = convert(Array{Float64}, genoX[:,2:end])  #first colum is ID
    X2           = convert(Array{Float64}, genoTrain2[:,2:end])  #first colum is ID
    genoX = 0 #release memory
    genoTrain2 = 0
    println("X is this size", size(X1),size(X2))
    Y           = convert(Array{Float64}, phenoTrain)
    println("Y is this size", size(Y))
    nTraits, nRecords , nMarkers   = size(Y,2), size(Y,1), size(X1,2)
    fileControl(nTraits,fixedRegSize)
    p1           = mean(X1,dims=1)./2.0
    sum2pq1      = sum(2*(1 .- p1).*p1)
    
    p2           = mean(X2,dims=1)./2.0
    sum2pq2      = sum(2*(1 .- p2).*p2)

    sum2pq       = sqrt.([sum2pq1; sum2pq2]*[sum2pq1; sum2pq2]')
    println(sum2pq)
    
    #priors
const    dfβ         = dfEffect + nTraits
const    scaleRes1    = varResidual1*(dfRes-2.0)/dfRes    
const    scaleRes2    = varResidual2*(dfRes-2.0)/dfRes    


    if varGenotypic==0.0
        covBeta  = fill([0.003 0;0 0.003],nRegions)
        Vb       = covBeta[1]
        else
        covBeta  = fill(varGenotypic./sum2pq,nRegions)
        Vb       = covBeta[1].*(dfβ-nTraits-1)
    end

    νS_e1           = scaleRes1*dfRes
    df_e            = dfRes
    νS_e2           = scaleRes2*dfRes

    
    #initial Beta values as "0"
    tempBetaMat     = zeros(Float64,nTraits,nMarkers)
    μ               = mean(Y,dims=1)    
    X1             .-= ones(Float64,nRecords)*2*p1    
    X2             .-= ones(Float64,nRecords)*2*p2
    
    MpM = []
    for j in 1:nMarkers
        tempM = Array{Float64}(nRecords,0)
        tempM = convert(Array{Float64},hcat(tempM,X1[:,j],X2[:,j]))
        this = tempM'tempM
        this[1,2]=this[2,1]=0.0
        MpM = push!(MpM,this)
    end
    nowM  = 0
    tempM = 0
    
    
    Ycorr = Y .- μ
            
    for iter in 1:chainLength
        #sample residual var
        R1 = sampleVarE(νS_e1,Ycorr[:,1],df_e,nRecords)
        R2 = sampleVarE(νS_e2,Ycorr[:,2],df_e,nRecords)
        Rmat = [R1 0;0 R2]
#        Ri = kron(inv(Rmat),eye(nRecords))
        Ri = inv(Rmat)

        Ycorr .+= μ
        for t in 1:nTraits
            rhs = sum(view(Ycorr,:,t))
            invLhs = 1.0/nRecords
            mean = rhs*invLhs
            μ[t] = rand(Normal(mean,sqrt(invLhs*Rmat[t,t])))
        end
        Ycorr .-= μ
        
        for r in 1:nRegions
            theseLoci = SNPgroups[r]
            regionSize = length(theseLoci)
            invB = inv(covBeta[r])
            for locus in theseLoci::UnitRange{Int64}
                sampleBeta_xiolei!(tempBetaMat,nTraits,X1,X2,Ri,locus,MpM,Ycorr,invB)
            end
            covBeta[r] = sampleCovBeta(dfβ,regionSize,Vb,tempBetaMat, theseLoci)
        end
        outputControl_xiolei(sum2pq,onScreen,iter,these2Keep,tempBetaMat,μ,covBeta,Rmat,fixedRegSize,nRegions)
    end
    GC.gc()
end

function sampleBeta_xiolei!(tempBetaMat,nTraits,X1,X2,Ri,locus,xpx,Ycorr,invB)
    Ycorr .+= [view(X1,:,locus).*view(tempBetaMat,1,locus) view(X2,:,locus).*view(tempBetaMat,2,locus)]
    rhs     = [view(X1,:,locus)'*view(Ycorr,:,1)*Ri[1];view(X2,:,locus)'*view(Ycorr,:,2)*Ri[4]]
    invLhs  = inv(xpx[locus].*Ri .+ invB)    
    meanBeta = invLhs*rhs
    tempBetaMat[:,locus] = rand(MvNormal(meanBeta,convert(Array,Symmetric(invLhs))))
    Ycorr .-= [view(X1,:,locus).*view(tempBetaMat,1,locus) view(X2,:,locus).*view(tempBetaMat,2,locus)]
end

function outputControl_xiolei(sum2pq,onScreen,iter,these2Keep,tempBetaMat,μ,covBeta,varE,fixedRegSize,nRegions)
    if iter in these2Keep
        out0 = open(pwd()*"/muOut$fixedRegSize", "a")
        writecsv(out0, μ)
        close(out0)
        for t in 1:2
            out1 = open(pwd()*"/beta"*"$t"*"Out$fixedRegSize", "a")
            writecsv(out1, tempBetaMat[t,:]')
            close(out1)
        end
        outCov = open(pwd()*"/covBetaOut$fixedRegSize", "a")
        printThis = [vcat(covBeta[r]...) for r in 1:nRegions]'
        writecsv(outCov, printThis)
        close(outCov)
        out3 = open(pwd()*"/varEOut$fixedRegSize", "a")
        writecsv(out3, varE)
        close(out3)
        if onScreen==true
            coVarBeta = cov(tempBetaMat')
            corBeta   = cor(tempBetaMat') 
            genCov    = sum2pq.*coVarBeta
            println("iter $iter \n coVarBeta (Overall): $coVarBeta \n genCov: $genCov \n corBeta: $corBeta \n varE: $varE \n")
        elseif onScreen==false
             @printf("iter %s\n", iter)
        end
    end
end
