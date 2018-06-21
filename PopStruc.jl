mutable struct Human
    strains_matrix::Array{Int64,2}
    Vector_time::Array{Int64,1}
    NumberStrains::Int64
    latenttime::Int64
    #=health::HEALTH
    swap::HEALTH #do we need  this? We can do a sequential atualization
    timeinstate::Int64
    statetime::Int64
    vaccinationStatus::Int64
    vaccineEfficacy::Float64
    WhoInf::Int64
    WentTo::HEALTH
    age::Int64
    group::Int64
    contact_group::Int64
    daily_contacts::Int64
    index::Int64
    Coverage::Float64
    NumberFails::Int64=#
   # Human() = new(SUSC,UNDEF,0,999,0,0.0,-1,UNDEF,-1,-1,-1,0,-1,0.0,0)
end

function setup_human(h)
    for i = 1:length(h)
        h[i] = create_human()
    end
end

function create_human()
    h = Human(zeros(Int64,P.matrix_strain_lines,P.sequence_size),zeros(Int64,P.matrix_strain_lines),0)
    #h.strains_matrix = zeros(Int64,P.matrix_strain_lines,P.sequence_size)
    return h
end


function setup_demographic(h::Array{Human},P::InfluenzaParameters)

    dist,AgeMin,AgeMax = distribution_age()

    for i = 1:17
        
        for aux2 = 1:2
            aux = 0
            while aux == 0
                j = rand(1:P.grid_size_human)
                if h[j].age < 0 
                    h[j].age = rand(AgeMin[i] : AgeMax[i])
                    if i <= 15
                        h[j].contact_group = i
                    else h[j].contact_group = 15
                    end
                    aux = 1
                end
            end
        end

    end 

    for i = 1:P.grid_size_human
        h[i].index = i
        if h[i].age < 0 
            rn = rand()
            g = findfirst(x -> rn <= x, dist)
            h[i].age = rand(AgeMin[g] : AgeMax[g])

            if g <= 15
                h[i].contact_group = g
            else h[i].contact_group = 15
            end
        end

        ###Group 1 - young child, 2 - school child, 3 - working adult, 4 - elderly
        if h[i].age<=4
            h[i].group = 1
        else
            if h[i].age<=19
            h[i].group = 2
            else
                if h[i].age<= 59
                    h[i].group = 3
                else h[i].group = 4
                end
            end
        end

        h[i].Coverage = settingCov(h[i])

    end

end

function settingCov(h::Human)

    if h.age<=4
        Coverage = 0.26
    else
        if h.age<= 49
            Coverage = 0.23
        else
            if h.age<= 64
                Coverage = 0.38
            else Coverage = 0.70
            end
        end
    end
    return Coverage
end


function ContactMatrixFunc()

    ContactMatrix = [0.2287227639 0.0578904992 0.0147540984 0.0100089736 0.0190779014 0.037582766 0.0554200656 0.0437937938 0.0462266543 0.0236936358 0.0184025802 0.0231225914 0.0311685233 0.0212494688 0.0201779473;
    0.3440471891 0.4866344605 0.0763661202 0.0301649755 0.0345345345 0.0714352327 0.1177361029 0.1086086086 0.093551931 0.0491321517 0.0539745779 0.0524945318 0.0557004029 0.0569485763 0.0478233238;
    0.3872637535 0.5714170692 0.5845628415 0.105957065 0.0548489666 0.0913923342 0.1621394332 0.1587420754 0.1551593003 0.0899990816 0.1026370708 0.0834288095 0.0850912539 0.0820229494 0.0986653956;
    0.4126640183 0.5926731079 0.6730874317 0.6061986609 0.1655184596 0.1379278187 0.1983853334 0.1941107774 0.2220907631 0.1799981633 0.1611648644 0.1249869805 0.1146006163 0.1053973651 0.1525262154;
    0.4493800409 0.6101449275 0.6949453552 0.6866155864 0.4279279279 0.2639186795 0.2519552603 0.241324658 0.291219471 0.2458444302 0.2378106621 0.1792521612 0.1570277317 0.1399631676 0.1839847474;
    0.5237751294 0.652173913 0.7161202186 0.7236142749 0.5620915033 0.4486617551 0.3590110167 0.3191524858 0.3590805375 0.322159978 0.3357996585 0.2558066868 0.2245792842 0.2059781839 0.2353034636;
    0.6348862405 0.7144122383 0.7481557377 0.7523296749 0.65527292 0.5724144363 0.5182070473 0.4225058392 0.4519563931 0.4088529709 0.4162398027 0.3453806895 0.3200995497 0.2838929027 0.2902764538;
    0.7340796918 0.791626409 0.8023907104 0.7966452682 0.722045575 0.6619416208 0.6553696073 0.5838338338 0.5565790586 0.5072091101 0.4973439575 0.4229767733 0.4124200047 0.3817821221 0.3541468065;
    0.7954736969 0.8639291465 0.8740437158 0.8515910817 0.7824589295 0.736174578 0.7554452948 0.7181348015 0.7006676244 0.6194324548 0.5895465756 0.5163003854 0.5003555345 0.4792463522 0.4401016841;
    0.8345973276 0.8972624799 0.9204234973 0.9178573894 0.8578872991 0.8051851161 0.8231435539 0.8019686353 0.805121271 0.756359629 0.6913299184 0.5997291949 0.5699217824 0.5429947585 0.5290753098;
    0.8760081859 0.9234299517 0.9460382514 0.9519569269 0.9194488606 0.8811899655 0.8819275082 0.8677844511 0.8823628835 0.8526035449 0.8022196927 0.7022185189 0.6420952832 0.6102847429 0.5967588179;
    0.9141687733 0.9438808374 0.9607240437 0.9704562711 0.9528351881 0.9308962044 0.9244806997 0.9051551552 0.9129552945 0.9038479199 0.8839878581 0.8314758879 0.743422612 0.6843745573 0.6561804894;
    0.9461899603 0.9631239936 0.9705601093 0.979153724 0.9685567921 0.959153222 0.9572786141 0.9432766099 0.9416884983 0.9339700615 0.9279074179 0.9002187272 0.8497274236 0.7788638617 0.7316491897;
    0.9695437583 0.9768921095 0.9801912568 0.9861945192 0.9777424483 0.9732351021 0.9746867379 0.9708041375 0.9620552692 0.9496739829 0.9462151394 0.9396937819 0.9156198151 0.8810029749 0.8156974897;
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ]
    
    return ContactMatrix
    
end 
