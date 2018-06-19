function mutation(Original_Strain1::Array{Int64,1},P::InfluenzaParameters,t::Int64,CumProb::Array{Float64,1},n::Int64)
   
    j = n-1 ###strain index
    Matrix_Of_Strains = zeros(Int64,P.matrix_strain_lines,P.sequence_size) ###Returning matrix with strains
    Matrix_Of_Strains[1,:] = Original_Strain1 #setting the first one as the original infection
    Time_Strain = zeros(Int64,P.matrix_strain_lines) ### the time in which the strain was generated
    Time_Strain[1] = 1 ### the original strain is generated at time 1
    

    while j != n ##I want to do this until there is no more strains to search
        j+=1 ##going to the next strain
        Time_Vector_Size = zeros(Int64,P.max_infectious_period) ##This vector saves the number of sites changing at that time step
        Time_Matrix = zeros(Int64,100,P.max_infectious_period) ##This matrix saves which sites changed at a given time step
        for i = 1:P.sequence_size ###running the sites
            if Time_Strain[j]<t    
                if rand() < CumProb[t-Time_Strain[j]] ##checking if it will change
                    time = rand(Time_Strain[j]+1:t) ###selecting a random time step to change it
                    Time_Vector_Size[time]+=1 ###increasing the number of sites that changed at time step "time"
                    Time_Matrix[Time_Vector_Size[time],time] = i ##saving the site
                end
            end
        end

        for time = Time_Strain[j]:t #running the time (each strain is able to create t-Time_Strain[j] strains)
            if Time_Vector_Size[time] > 0 ##if it has sites that changed
                n+=1
                Matrix_Of_Strains[n,:] = Matrix_Of_Strains[j,:] ###Copying the generating strain
                for time_count = 1:Time_Vector_Size[time] #running the sites that changed

                    ##this loop guarantees that the new state is different from the old one
                    change::Int64 = Matrix_Of_Strains[n,Time_Matrix[time_count,time]] 
                    while change == Matrix_Of_Strains[n,Time_Matrix[time_count,time]]
                        change = rand(1:20)
                    end
                    ###########################
                    Matrix_Of_Strains[n,Time_Matrix[time_count,time]] = change ##changing the site's state
                end
                Time_Strain[n] = time #saving the time the strain appeared
            end
        end ###Close for time
    end

    return Matrix_Of_Strains,Time_Strain,n ##return the strain matrix, when they were generated and how many they are

end


function Calculating_Distance(humans::Array{Human},P::InfluenzaParameters)
    n::Int64 = 0
    for i = 1:P.grid_size_human
        n+=humans[i].NumberStrains
    end

    A = zeros(Int64,P.sequence_size,n)
    count::Int64 = 0
    for i = 1:P.grid_size_human
        for j = 1:humans[i].NumberStrains
            count+=1
            A[:,count] = humans[i].strains_matrix[j,:]
        end
    end

    DistMatrix = zeros(Int64,n,n)
    for i = 1:(n-1)
        for j = i+1:n
           DistMatrix[i,j] = Calculating_Distance_Two_Strains(A[:,i],A[:,j])
        end
    end

    return DistMatrix
end

function Calculating_Distance_Two_Strains(A::Array{Int64,1},B::Array{Int64,1})
    soma::Int64 = 0
    for  k = 1:length(A)
        if A[k] != B[k]
            soma+=1
        end
    end
    return soma
end