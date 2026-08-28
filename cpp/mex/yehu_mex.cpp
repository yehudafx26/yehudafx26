#include "mex.hpp"
#include "mexAdapter.hpp"
#include <vector>
#include <string>
#include "yehu.h"
#include <chrono>
#include <iostream>
#include <iomanip>

/**
 * MEX WRAPPER FUNCTION FOR YEHU SIMULATIONS
 * Champ C. Darabundit and Zhen Zhang, CAML 2026
 */
class MexFunction : public matlab::mex::Function {
    public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
        {
            checkArguments(outputs, inputs);
            
            // ── Nf x 1 column vectors ──────────────────────────────────────────
            std::vector<double> vb = toVector(inputs[0]);
            std::vector<double> fn = toVector(inputs[1]);
            std::vector<double> xf = toVector(inputs[2]);
            std::vector<double> ff = toVector(inputs[3]);
            
            const size_t Nf = vb.size();
            
            // ── 1 x 8 row vectors ─────────────────────────────────────────────
            std::vector<double> fStringParams = toVector(inputs[4]);
            std::vector<double> cStringParams = toVector(inputs[5]);

            // ── bodyParams M x 3 → Mq, Kq, Rq ─────────────────────────────
            matlab::data::TypedArray<double> bodyArr = inputs[6];
            const auto bodyDims = bodyArr.getDimensions();
            const size_t M = bodyDims[0];
            
            std::vector<double> Mq(M), Kq(M), Rq(M);
            for (size_t i = 0; i < M; ++i) {
                Mq[i] = bodyArr[i][0];
                Kq[i] = bodyArr[i][1];
                Rq[i] = bodyArr[i][2];
            }
            
            // ── 1 x 5 row vector ──────────────────────────────────────────────
            std::vector<double> fingerParams = toVector(inputs[7]);
            
            // ── 1 x 4 row vector ──────────────────────────────────────────────
            std::vector<double> bowParams = toVector(inputs[8]);
            
            // ── 1 x 6 row vector ──────────────────────────────────────────────
            std::vector<double> frictionParams = toVector(inputs[9]);
            
            // ── SR scalar ────────────────────────────────────────────────────
            matlab::data::TypedArray<double> srArr = inputs[10];
            const double SR = srArr[0];
            
            // // PARAMETER INITIALIZATION
            sim.setString1Params( fStringParams[0], fStringParams[1], fStringParams[2], fStringParams[3], fStringParams[4], fStringParams[5], fStringParams[6], fStringParams[7] );
            sim.setString2Params( cStringParams[0], cStringParams[1], cStringParams[2], cStringParams[3], cStringParams[4], cStringParams[5], cStringParams[6], cStringParams[7] );
            sim.setBodyParameters(Mq, Kq, Rq);
            sim.setFingerParams( fingerParams[0],fingerParams[1], fingerParams[2], fingerParams[3], fingerParams[4] ) ;
            sim.setBowParams( bowParams[0],bowParams[1], bowParams[2], bowParams[3] );
            sim.setFrictionParams( frictionParams[0], frictionParams[1], frictionParams[2], frictionParams[3], frictionParams[4], frictionParams[5] );
            sim.prepareToPlay( SR );
            
            // COMPUTE HERE
            std::vector<double> out(Nf, 0.0);
            std::vector<double> fbri(Nf, 0.0);
            auto tic = std::chrono::high_resolution_clock::now();
            sim.process(out, fbri, vb, fn, xf, ff);
            std::chrono::duration<double> toc = std::chrono::high_resolution_clock::now() - tic;
            std::cout << "Elpased time is " << std::setprecision(8) << toc.count() << " seconds." << std::endl;
            
            
            // ── Output: Nf x 1 column vector ──────────────────────────────────
            matlab::data::TypedArray<double> outArr =
            factory.createArray<double>({ Nf, 1 }, out.data(), out.data() + Nf);
             matlab::data::TypedArray<double> fbriArr =
            factory.createArray<double>({ Nf, 1 }, fbri.data(), fbri.data() + Nf);
            outputs[0] = outArr;
            outputs[1] = fbriArr;
        }
        
        void checkArguments(matlab::mex::ArgumentList outputs,
        matlab::mex::ArgumentList inputs)
        {
            if (outputs.size() != 2) {
                throwError("yehu:badNumOutputs", "Expected 2 output arguments.");
            }
            
            if (inputs.size() != 11) {
                throwError("yehu:badNumInputs",
                    "Expected 11 input arguments: "
                    "vb, fn, xf, ff, fStringParams, cStringParams, "
                    "fingerParams, bowParams, frictionParams, bodyParams, SR.");
            }
                
            // vb, fn, xf, ff — Nf x 1, all the same length
            const size_t Nf = inputs[0].getDimensions()[0];
            const char* vecNames[4] = { "vb", "fn", "xf", "ff" };
            
            for (int i = 0; i < 4; ++i) {
                const auto dims = inputs[i].getDimensions();
                if (dims.size() != 2 || dims[1] != 1) {
                    throwError("yehu:badInput", std::string(vecNames[i]) + " must be an Nf x 1 column vector.");
                }
                if (dims[0] != Nf) {
                    throwError("yehu:badInput", std::string(vecNames[i]) +" must have the same length as vb (Nf x 1).");
                }
            }
                
                // fStringParams — 1 x 8
                checkSize(inputs[4], 8, "fStringParams");
                
                // cStringParams — 1 x 8
                checkSize(inputs[5], 8, "cStringParams");

                // bodyParams — M x 3
                const auto bodyParamDims = inputs[6].getDimensions();
                if (bodyParamDims.size() != 2 || bodyParamDims[1] != 3) {
                    throwError("yehu:badInput", "bodyParams must be an M x 3 matrix.");
                }

                // fingerParams — 1 x 5
                checkSize(inputs[7], 5, "fingerParams");
                
                // bowParams — 1 x 4
                checkSize(inputs[8], 4, "bowParams");
                
                // frictionParams — 1 x 8
                checkSize(inputs[9], 6, "frictionParams");
                
                // SR — scalar double
                const auto srDims = inputs[10].getDimensions();
                if (srDims[0] != 1 || srDims[1] != 1) {
                    throwError("yehu:badInput", "SR must be a scalar double.");
                }
        }

        private:
        yehu::yehu<double> sim;
        matlab::data::ArrayFactory factory;
        std::shared_ptr<matlab::engine::MATLABEngine> engine = getEngine();

        // ── Argument validation ───────────────────────────────────────────────
        
            
        void checkSize(const matlab::data::Array& arr, size_t expectedElements, const std::string& name)
        {
            if (arr.getNumberOfElements() != expectedElements) {
                throwError("yehu:badInput",
                    name + " must have " +
                    std::to_string(expectedElements) + " elements.");
                }
        }
                
        std::vector<double> toVector(const matlab::data::Array& arr) {
            const matlab::data::TypedArray<double> typed = arr;
            return std::vector<double>(typed.begin(), typed.end());
        }
                
        void throwError(const std::string& id, const std::string& msg) {
            engine->feval(u"error", 0, std::vector<matlab::data::Array>({ factory.createScalar(id), factory.createScalar(msg)}));
        }

};