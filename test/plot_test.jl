using AdiabaticEvolution
using Plots
using Test

@testset "Plotting Tests" begin
    # Generate some test data
    times = range(0, 10, length=100)
    energies = sin.(times) .+ 2.0  # Simple sinusoidal energy evolution
    
    # Test that the plotting function returns a plot object
    p = plot_adiabatic_evolution(times, energies)
    @test p isa Plots.Plot
    
    # Test that the plot has the correct attributes
    @test p.attr[:xaxis][:label] == "Time"
    @test p.attr[:yaxis][:label] == "Energy"
    @test p.attr[:title] == "Adiabatic Evolution"
    
    # Save the plot for visual inspection
    savefig(p, "test_plot.png")
    @test isfile("test_plot.png")
    
    # Clean up
    rm("test_plot.png")
end 