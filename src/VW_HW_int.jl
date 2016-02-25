module VictorWarhem_solutions_PBS1

	using PyPlot
	using FastGaussQuadrature
	using Roots
 	using Sobol
 	using PyPlot
 	using Distributions
 	using Calculus

 	function question_1a(n)

 		# We define the space and the demand function:

 		p=linspace(0,4,n)
		q(p)=2*p .^(-0.5)

		# We plot:

		plot(p,q(p),color="red")
		q1=linspace(q(1),q(1),n)
		q2=linspace(q(4),q(4),n)
		plot(p,q1,color="blue")
		plot(p,q2,color="blue")
		title("Demand function and optimal quantities")
		xlabel("p")
		ylabel("q(p)")
	end

	question_1a(10)
	question_1a(100)
	question_1a(1000)

	function question_1b(n)
    
        p=linspace(0,4,n)
		q(p)=2*p .^(-0.5)

		# We create a two columns matrix containing the nodes and weights calculated be the gauss legendre method:
    
		nodes,weights= gausslegendre(n)

		# We define a function y:

		y=[]
		for i in 1:n
		push!(y,weights[i]*q(1.5*nodes[i]+2.5))
		end

		# We finish the approximation:

		approxgausslegendre=1.5*sum(y,1)

		# We plot:

		
		plot(p,q(1.5*p+2.5),color="blue")
		title("Change in Consumer Surplus with Gauss-Legendre")
		xlabel("nodes")
		ylabel("q(1.5*nodes+2.5)")

		println("The numerical change in Consumer Surplus with the Gauss Legendre integration method with n = $n is $(round(approxgausslegendre,5))")

		println("")
	end

	question_1b(10)
	question_1b(100)
	question_1b(1000)

	# We use the Monte Carlo integration method:

	function question_1c(n)
    
        p=linspace(0,4,n)
		q(p)=2*p^(-0.5)
    
        nodes=3rand()+1
    
		G=[]
		for i in 1:n
		push!(G,q(nodes))
		end

		approxmontecarlo=3sum(G,1)/n

		qu(p)=2*p .^(-0.5)
		nodes=3*rand(n)+1
		plot(p,qu(p),color="blue")
		scatter(nodes, qu(nodes), color ="black" )
		title("Change in Consumer Surplus with Monte Carlo")
		xlabel("nodes")
		ylabel("qu(nodes)")

		println("The numerical change in Consumer Surplus with the Monte Carlo integration method with n = $n is $(round(approxmontecarlo,5))")

		println("")
	end

	question_1c(10)
	question_1c(100)
	question_1c(1000)

	# We replace the random nodes by a Sobol sequence:

	using Sobol

	function question_1d(n)
        
        s=SobolSeq(1,1,4)

		p=linspace(0,4,n)
		q(p)=2*p^(-0.5)

		H=[]
		for i in 1:n
		nodes=hcat(next(s))
		push!(H,q(nodes))
		end

		approxquasimontecarlo=3sum(H,1)/n

		qu(p)=2*p .^(-0.5)
		nodes=3*rand(n)+1
		plot(p,qu(p),color="blue")
		scatter(nodes, qu(nodes), color ="black" )
		title("Change in Consumer Surplus with Quasi Monte Carlo")
		xlabel("nodes")
		ylabel("qu(nodes)")

		println("The numerical change in Consumer Surplus with the Quasi Monte Carlo integration method with n = $n is $(round(approxquasimontecarlo,5))")

		println("")
	end

	question_1d(10)
	question_1d(100)
	question_1d(1000)

	function question_2a(n)

		Sigma=[0.02 0.01; 0.01 0.01]
		Halfsigma=chol(Sigma,Val{:L})

		# First we use the Gauss-Hermite Integration method
		
		nodes1, weights1 = gausshermite(n)
		nodes2, weights2 = gausshermite(n)

		nodes = hcat(kron(ones(n),nodes1),kron(nodes2,ones(n)))

		# we define the weights and what's inside the function after changing the variable:

		w = kron(weights1,weights2) / sqrt(pi)
		g =  nodes * Halfsigma + zeros(n*n,2)
		
		# We calculate p such that nodes always respect the market clearing condition

		p = []

		for i = 1:n*n
		function dd(p)
				exp(g[i,1])*p^(-1) .+ exp(g[i,2]) .* p^(-0.5) - 2
		end
		push!(p,fzeros(dd, [0,3]))
		end

		# We find the expectation
 
    	resultexp=sum(w .*p)
    
    	println("I find as expectancy with n = $n $(round(resultexp,5)).")

    	println("No results for the variance.")

    	println("")
	end

	question_2a(10)

 	function question_2b(n)

 		Sigma=[0.02 0.01; 0.01 0.01]
		Halfsigma=chol(Sigma,Val{:L})

  		# Now we use the lognormal distribution to generate the nodes and the weights by assuming a log function that is normal implies a function that is lognormal.
 
 		log1 = LogNormal(0.0,0.02)
 		log2 = LogNormal(0.0,0.01)
 		nodes1 = rand(log1,n)
 		nodes2 = rand(log2,n)
 
 		nodes = hcat(kron(ones(n),nodes1),kron(nodes2,ones(n)))

		# we define the weights and what's inside the function after changing the variable:

		w = ones(n*n) ./ n
		g =  nodes * Halfsigma + zeros(n*n,2)
		
		# We calculate p such that nodes always respect the market clearing condition

		p = []

		for i = 1:n*n
		function dd(p)
				exp(g[i,1])*p^(-1) .+ exp(g[i,2]) .* p^(-0.5) - 2
		end
		push!(p,fzeros(dd, [0,3]))
		end

		# We find the expectation
 
    	resultexp=sum(w .*p)
    
    	println("I find as expectancy with n = $n $(round(resultexp,5)), which is not a good approximation, as everyone can tell.")

    	println("No results for the variance.")

    	println("")
 	end
 
 question_2b(10)

 	function runall(n=10)
 		println("running all questions of HW-integration:")
 		println("results of question 1:")
 		question_1a(n)
 		question_1b(n)
 		question_1c(n)
 		question_1d(n)
 		println("results of question 2:")
 		question_2a(n)
 		question_2b(n)
 		println("end of HW-integration")
 		println("")
 	end
 
end



