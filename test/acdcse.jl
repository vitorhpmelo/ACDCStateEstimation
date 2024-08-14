# Test problem functions

# The purpose of the tests contained in this file is to detect if anything has accidentally
# changed in the problem functions. Accordingly, only termination status and objective value
# are tested.
# For testing specific features, it is better to write ad-hoc tests in separate files.

@testset "ACDCSE" begin

    data = _ACDCSE.quickget_case5()  
    pf_result = _PMMCDC.solve_mcdcopf(data, _PM.ACPPowerModel, nlp_optimizer)
    _ACDCSE.get_dc_power!(pf_result, data)

    @testset "ACDCSE ACP - no measurement conversion needed" begin 
        
        @testset "case5_2grids_MC - noiseless SE" begin
            add_noise = false
            _ACDCSE.powerflow2measurements!(data, pf_result, σ_dict, sample_error = add_noise, 
                    measurements = ["vm", "va", "p_fr", "q_fr", "p_to", "q_to", "vdcm", "pg", "pd", "qg", "qd", "i_dcgrid_fr", "i_dcgrid_to"])

            _ACDCSE.prepare_data_for_se_default!(data, exceptions = [1]) 
            se_res = _ACDCSE.solve_acdcse(data, _PM.ACPPowerModel, nlp_optimizer)

            @test se_res["termination_status"] == _PMMCDC.LOCALLY_SOLVED
            @test isapprox(se_res["objective"], 0.0; atol=1e-4)
        end

        @testset "case5_2grids_MC - noisy SE" begin
            add_noise = true
            _ACDCSE.powerflow2measurements!(data, pf_result, σ_dict, sample_error = add_noise, 
                    measurements = ["vm", "va", "p_fr", "q_fr", "p_to", "q_to", "vdcm", "pg", "pd", "qg", "qd", "i_dcgrid_fr", "i_dcgrid_to"])

            _ACDCSE.prepare_data_for_se_default!(data, exceptions = [1]) 
            se_res = _ACDCSE.solve_acdcse(data, _PM.ACPPowerModel, nlp_optimizer)

            @test se_res["termination_status"] == _PMMCDC.LOCALLY_SOLVED
            @test isapprox(se_res["objective"], 83.534; atol=1e-3)
        end

        @testset "ACDCSE ACP - with measurement conversions" begin 
        
            @testset "case5_2grids_MC - noiseless SE" begin
                add_noise = false
                _ACDCSE.powerflow2measurements!(data, pf_result, σ_dict, sample_error = add_noise, 
                        measurements = ["vm", "pg", "pd", "qg", "qd", "vdcm", "p_dc_fr", "p_dc_to", "i_dcgrid_fr", "i_dcgrid_to"])
                        # TODO !!! ↑ add here some current magnitudes!
                _ACDCSE.prepare_data_for_se_default!(data, exceptions = [1]) 
                se_res = _ACDCSE.solve_acdcse(data, _PM.ACPPowerModel, nlp_optimizer)
    
                @test se_res["termination_status"] == _PMMCDC.LOCALLY_SOLVED
                @test isapprox(se_res["objective"], 0.0; atol=1e-4)
            end
    
            # @testset "case5_2grids_MC - noisy SE" begin
            #     add_noise = true
            #     _ACDCSE.powerflow2measurements!(data, pf_result, σ_dict, sample_error = add_noise, 
            #             measurements = ["vm", "va", "p_fr", "q_fr", "p_to", "q_to", "vdcm", "pg", "pd", "qg", "qd", "i_dcgrid_fr", "i_dcgrid_to"])
    
            #     _ACDCSE.prepare_data_for_se_default!(data, exceptions = [1]) 
            #     se_res = _ACDCSE.solve_acdcse(data, _PM.ACPPowerModel, nlp_optimizer)
    
            #     @test se_res["termination_status"] == _PMMCDC.LOCALLY_SOLVED
            #     display(se_res["objective"])
            #     # @test isapprox(se_res["objective"], 0.0; atol=1e-4)
            # end
        end
    end

end