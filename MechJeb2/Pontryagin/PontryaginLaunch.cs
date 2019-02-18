using System;
using System.Collections.Generic;
using System.Threading;
using UnityEngine;

namespace MuMech {
    public class PontryaginLaunch : PontryaginBase {
        public PontryaginLaunch(MechJebCore core, double mu, Vector3d r0, Vector3d v0, Vector3d pv0, Vector3d pr0, double dV) : base(core, mu, r0, v0, pv0, pr0, dV)
        {
        }

        public bool omitCoast;

        double rTm;
        double vTm;
        double gammaT;
        double incT;
        double smaT;
        double eccT;
        double LANT;
        double ArgPT;
        Vector3d hT;

        // 5-constraint PEG with fixed LAN
        public void flightangle5constraint(double rTm, double vTm, double gamma, Vector3d hT)
        {
            QuaternionD rot = Quaternion.Inverse(Planetarium.fetch.rotation);
            this.rTm = rTm / r_scale;
            this.vTm = vTm / v_scale;
            this.gammaT = gamma;
            this.hT = rot * hT / r_scale / v_scale;
            bcfun = flightangle5constraint;
        }

        private void flightangle5constraint(double[] yT, double[] z, bool terminal)
        {
            Vector3d rf = new Vector3d(yT[0], yT[1], yT[2]);
            Vector3d vf = new Vector3d(yT[3], yT[4], yT[5]);
            Vector3d pvf = new Vector3d(yT[6], yT[7], yT[8]);
            Vector3d prf = new Vector3d(yT[9], yT[10], yT[11]);
            Vector3d hf = Vector3d.Cross(rf, vf);

            Vector3d hmiss = hf - hT;

            // 5 constraints
            if (!terminal)
            {
                z[0] = ( rf.magnitude * rf.magnitude - rTm * rTm ) / 2.0;
                z[1] = Vector3d.Dot(rf, vf) - rf.magnitude * vf.magnitude * Math.Sin(gammaT);
                z[2] = hmiss[0];
                z[3] = hmiss[1];
                z[4] = hmiss[2];

                // transversality - free argp
                z[5] = Vector3d.Dot(Vector3d.Cross(prf, rf) + Vector3d.Cross(pvf, vf), hT);
            }
            else
            {
                z[0] = hmiss.magnitude;
                z[1] = z[2] = z[3] = z[4] = z[5] = 0.0;
            }
        }

        // 4-constraint PEG with free LAN
        public void flightangle4constraint(double rTm, double vTm, double gamma, double inc)
        {
            //Debug.Log("call stack: + " + Environment.StackTrace);
            this.rTm = rTm / r_scale;
            this.vTm = vTm / v_scale;
            //Debug.Log("4constraint vTm = " + vTm + " v_scale = " + v_scale + " vTm_bar = " + this.vTm );
            //Debug.Log("4constraint rTm = " + rTm + " r_scale = " + r_scale + " rTm_bar = " + this.rTm );
            this.gammaT = gamma;
            this.incT = inc;
            bcfun = flightangle4constraint;
        }

        private void flightangle4constraint(double[] yT, double[] z, bool terminal)
        {
            Vector3d rf = new Vector3d(yT[0], yT[1], yT[2]);
            Vector3d vf = new Vector3d(yT[3], yT[4], yT[5]);
            Vector3d pvf = new Vector3d(yT[6], yT[7], yT[8]);
            Vector3d prf = new Vector3d(yT[9], yT[10], yT[11]);

            Vector3d n = new Vector3d(0, -1, 0);  /* angular momentum vectors point south in KSP and we're in xzy coords */
            Vector3d rn = Vector3d.Cross(rf, n);
            Vector3d vn = Vector3d.Cross(vf, n);
            Vector3d hf = Vector3d.Cross(rf, vf);

            if (!terminal)
            {
                z[0] = ( rf.magnitude * rf.magnitude - rTm * rTm ) / 2.0;
                z[1] = ( vf.magnitude * vf.magnitude - vTm * vTm ) / 2.0;
                z[2] = Vector3d.Dot(n, hf) - hf.magnitude * Math.Cos(incT);
                z[3] = Vector3d.Dot(rf, vf) - rf.magnitude * vf.magnitude * Math.Sin(gammaT);
                z[4] = rTm * rTm * ( Vector3d.Dot(vf, prf) - vTm * Math.Sin(gammaT) / rTm * Vector3d.Dot(rf, prf) ) -
                    vTm * vTm * ( Vector3d.Dot(rf, pvf) - rTm * Math.Sin(gammaT) / vTm * Vector3d.Dot(vf, pvf) );
                z[5] = Vector3d.Dot(hf, prf) * Vector3d.Dot(hf, rn) + Vector3d.Dot(hf, pvf) * Vector3d.Dot(hf, vn);
            }
            else
            {
                double hTm = rTm * vTm * Math.Cos(gammaT);

                z[0] = hf.magnitude - hTm;
                z[1] = z[2] = z[3] = z[4] = z[5] = 0.0;
            }
        }

        public void target3constraint(double sma, double ecc, double inc)
        {
            this.smaT = sma / r_scale;
            this.eccT = ecc;
            this.incT = inc;
            bcfun = target3constraint;
        }

        // has a singularity for e == 0 which is only fixable by going to flightangle4constraint
        // FIXME: singularity for i == 0, rot coordinates?
        private void target3constraint(double[] yT, double[] z, bool terminal)
        {
            Vector3d rf = new Vector3d(yT[0], yT[1], yT[2]);
            Vector3d vf = new Vector3d(yT[3], yT[4], yT[5]);
            Vector3d pvf = new Vector3d(yT[6], yT[7], yT[8]);
            Vector3d prf = new Vector3d(yT[9], yT[10], yT[11]);

            Vector3d hf = Vector3d.Cross(rf, vf);
            Vector3d n = new Vector3d(0, -1, 0);  /* angular momentum vectors point south in KSP and we're in xzy coords */

            double hTm = Math.Sqrt( smaT * ( 1 - eccT * eccT ) );

            double smaf = 1.0 / ( 2.0 / rf.magnitude - vf.sqrMagnitude );
            double eccf = Math.Sqrt(1.0 - hf.sqrMagnitude / smaf);

            if (!terminal)
            {
                z[0] = smaT * ( 1 - eccT ) - ( smaf * ( 1 - eccf ) ); // PeA constraint
                z[1] = vf.sqrMagnitude / 2.0 - 1.0 / rf.magnitude + 1.0 / ( 2.0 * smaT ); // E constraint
                z[2] = Vector3d.Dot(n, hf.normalized) - Math.Cos(incT); // inc constraint
                // transversality
                z[3] = Vector3d.Dot(Vector3d.Cross(prf, rf) + Vector3d.Cross(pvf, vf), hf);
                z[4] = Vector3d.Dot(Vector3d.Cross(prf, rf) + Vector3d.Cross(pvf, vf), n);
                z[5] = Vector3d.Dot(prf, vf) - Vector3d.Dot(pvf, rf) / ( rf.magnitude * rf.magnitude * rf.magnitude );
            }
            else
            {
                z[0] = hf.magnitude - hTm;
                z[1] = z[2] = z[3] = z[4] = z[5] = 0.0;
            }
        }

        public void target4constraintArgPfree(double sma, double ecc, double inc, double LAN)
        {
            this.smaT = sma / r_scale;
            this.eccT = ecc;
            this.incT = inc;
            this.LANT = LAN;
            bcfun = target4constraintArgPfree;
        }

        private void target4constraintArgPfree(double[] yT, double[] z, bool terminal)
        {
            Vector3d rf = new Vector3d(yT[0], yT[1], yT[2]);
            Vector3d vf = new Vector3d(yT[3], yT[4], yT[5]);
            Vector3d pvf = new Vector3d(yT[6], yT[7], yT[8]);
            Vector3d prf = new Vector3d(yT[9], yT[10], yT[11]);

            Vector3d hf = Vector3d.Cross(rf, vf);
            Vector3d n = new Vector3d(0, -1, 0);  /* angular momentum vectors point south in KSP and we're in xzy coords */

            Vector3d hT = new Vector3d( Math.Sin(LANT) * Math.Sin(incT), -Math.Cos(LANT) * Math.Sin(incT), Math.Cos(incT) ) * Math.Sqrt(smaT * (1 - eccT*eccT));
            hT = -hT.xzy; // left handed coordinate system
            Vector3d hmiss = hf - hT;

            double smaf = 1.0 / ( 2.0 / rf.magnitude - vf.sqrMagnitude );
            double eccf = Math.Sqrt(1.0 - hf.sqrMagnitude / smaf);

            if (!terminal)
            {
                z[0] = smaT * ( 1 - eccT ) - ( smaf * ( 1 - eccf ) ); // PeA constraint
                z[2] = hmiss[0];
                z[2] = hmiss[1];
                z[3] = hmiss[2];
                // transversality
                z[4] = Vector3d.Dot(Vector3d.Cross(prf, rf) + Vector3d.Cross(pvf, vf), hf);
                z[5] = Vector3d.Dot(prf, vf) - Vector3d.Dot(pvf, rf) / ( rf.magnitude * rf.magnitude * rf.magnitude );
            }
            else
            {
                z[0] = hmiss.magnitude;
                z[1] = z[2] = z[3] = z[4] = z[5] = 0.0;
            }
        }

        public void target4constraintLANfree(double sma, double ecc, double inc, double ArgP)
        {
            this.smaT = sma / r_scale;
            this.eccT = ecc;
            this.incT = inc;
            this.ArgPT = ArgP;
            bcfun = target4constraintLANfree;
        }

        private void target4constraintLANfree(double[] yT, double[] z, bool terminal)
        {
            Vector3d rf = new Vector3d(yT[0], yT[1], yT[2]);
            Vector3d vf = new Vector3d(yT[3], yT[4], yT[5]);
            Vector3d pvf = new Vector3d(yT[6], yT[7], yT[8]);
            Vector3d prf = new Vector3d(yT[9], yT[10], yT[11]);

            Vector3d hf = Vector3d.Cross(rf, vf);
            Vector3d n = new Vector3d(0, -1, 0); // angular momentum vectors point south in KSP and we're in xzy coords

            Vector3d eccf = Vector3d.Cross(vf, hf) - rf / rf.magnitude; // ecc vector
            double smaf = 1.0 / ( 2.0 / rf.magnitude - vf.sqrMagnitude );
            double hTm = Math.Sqrt( smaT * ( 1 - eccT * eccT ) );

            if (!terminal)
            {
                z[0] = smaT * ( 1 - eccT ) - ( smaf * ( 1 - eccf.magnitude ) ); // PeA constraint
                z[1] = vf.sqrMagnitude / 2.0 - 1.0 / rf.magnitude + 1.0 / ( 2.0 * smaT ); // E constraint
                z[2] = Vector3d.Dot(n, hf) - hf.magnitude * Math.Cos(incT);
                z[3] = Vector3d.Dot(eccf, Vector3d.Cross(n, hf)) / eccf.magnitude / hf.magnitude - Math.Sin(incT) * Math.Cos(ArgPT);
                z[4] = Vector3d.Dot(Vector3d.Cross(prf, rf) + Vector3d.Cross(pvf, vf), n);
                z[5] = Vector3d.Dot(prf, vf) - Vector3d.Dot(pvf, rf) / ( rf.magnitude * rf.magnitude * rf.magnitude );
            }
            else
            {
                z[0] = hf.magnitude - hTm;
                z[1] = z[2] = z[3] = z[4] = z[5] = 0.0;
            }
        }

        public void target5constraint(double sma, double ecc, double inc, double LAN, double ArgP)
        {
            this.smaT = sma / r_scale;
            this.eccT = ecc;
            this.incT = inc;
            this.LANT = LAN;
            this.ArgPT = ArgP;
            bcfun = target5constraint;
        }

        private void target5constraint(double[] yT, double[] z, bool terminal)
        {
            Vector3d rf = new Vector3d(yT[0], yT[1], yT[2]);
            Vector3d vf = new Vector3d(yT[3], yT[4], yT[5]);
            Vector3d pvf = new Vector3d(yT[6], yT[7], yT[8]);
            Vector3d prf = new Vector3d(yT[9], yT[10], yT[11]);

            Vector3d hT = new Vector3d( Math.Sin(LANT) * Math.Sin(incT), -Math.Cos(LANT) * Math.Sin(incT), Math.Cos(incT) ) * Math.Sqrt(smaT * (1 - eccT*eccT));
            hT = -hT.xzy; // left handed coordinate system
            // FIXME: Vector3d.cross(n, hf) is the node vector pointing in LAN dir, another ArgPT rot around hT would give eT direction
            Vector3d right = new Vector3d(1, 0, 0);
            Vector3d eT = Quaternion.Euler((float)LANT, (float)incT, (float)ArgPT) * right * (float) eccT;

            if (Math.Abs(hT[1]) <= 1e-6) // handle singularity
            {
                rf = rf.Reorder(231);
                vf = vf.Reorder(231);
                prf = prf.Reorder(231);
                pvf = pvf.Reorder(231);
                hT = hT.Reorder(231);
                eT = eT.Reorder(231);
            }

            Vector3d hf = Vector3d.Cross(rf, vf);
            Vector3d ef = - ( rf.normalized + Vector3d.Cross(hf, vf) );
            Vector3d hmiss = hf - hT;
            Vector3d emiss = ef - eT;
            double trans = Vector3d.Dot(prf, vf) - Vector3d.Dot(pvf, rf) / ( rf.magnitude * rf.magnitude * rf.magnitude );

            if (!terminal)
            {
                z[0] = hmiss[0];
                z[1] = hmiss[1];
                z[2] = hmiss[2];
                z[3] = emiss[0];
                z[4] = emiss[2];
                z[5] = trans;
            }
            else
            {
                z[0] = hmiss.magnitude;
                z[1] = z[2] = z[3] = z[4] = z[5] = 0.0;
            }
        }

        public override void Bootstrap(double t0)
        {
            // build arcs off of ksp stages, with coasts
            List<Arc> arcs = new List<Arc>();
            for(int i = 0; i < stages.Count; i++)
            {
                //if (i != 0)
                    //arcs.Add(new Arc(new Stage(this, m0: stages[i].m0, isp: 0, thrust: 0)));
                arcs.Add(new Arc(this, stages[i]));
            }

            arcs[arcs.Count-1].infinite = true;

            // allocate y0
            y0 = new double[arcIndex(arcs, arcs.Count)];

            // update initial position and guess for first arc
            double ve = g0 * stages[0].isp;
            tgo = ve * stages[0].startMass / stages[0].startThrust * ( 1 - Math.Exp(-dV/ve) );
            tgo_bar = tgo / t_scale;

            // initialize overall burn time
            y0[0] = tgo_bar;
            y0[1] = 0;

            UpdateY0(arcs);

            // seed continuity initial conditions
            yf = new double[arcs.Count*13];
            multipleIntegrate(y0, yf, arcs, initialize: true);

            //Debug.Log("optimizer hT = " + hT.magnitude * r_scale * v_scale);
            //Debug.Log("r_scale = " + r_scale);
            //Debug.Log("v_scale = " + v_scale);
            //Debug.Log("rTm = " + rTm * r_scale + " " + rTm);
            //Debug.Log("vTm = " + vTm * v_scale + " " + vTm);

            //for(int j = 0; j < y0.Length; j++)
            //    Debug.Log("bootstrap - y0[" + j + "] = " + y0[j]);

            //Debug.Log("running optimizer");

            if ( !runOptimizer(arcs) )
            {
                //for(int k = 0; k < y0.Length; k++)
                    //Debug.Log("failed - y0[" + k + "] = " + y0[k]);
                Fatal("Target is unreachable even with infinite ISP");
                return;
            }

            //Debug.Log("optimizer done");

            Solution new_sol = new Solution(t_scale, v_scale, r_scale, t0);
            multipleIntegrate(y0, new_sol, arcs, 10);

            for(int i = arcs.Count - 1; i >= 0; i--)
            {
                if ( new_sol.tgo(new_sol.t0, i) < 1 )
                {
                    /* burn is less than one second, delete it */
                    RemoveArc(arcs, i, new_sol);
                }
            }

            //for(int k = 0; k < y0.Length; k++)
            //   Debug.Log("y0[" + k + "] = " + y0[k]);

            //Debug.Log("arcs.Count = " + arcs.Count);
            //Debug.Log("segments.Count = " + new_sol.segments.Count);

            bool insertedCoast = false;

            if (arcs.Count > 1 && !omitCoast)
            {
                InsertCoast(arcs, arcs.Count-1, new_sol);
                insertedCoast = true;
            }
            arcs[arcs.Count-1].infinite = false;

            //Debug.Log("running optimizer3");

            if ( !runOptimizer(arcs) )
            {
                Fatal("Target is unreachable");
                //for(int k = 0; k < y0.Length; k++)
                    //Debug.Log("failed - y0[" + k + "] = " + y0[k]);
                //Debug.Log("optimizer failed6");
                y0 = null;
                return;
            }

            //Debug.Log("optimizer done");

            new_sol = new Solution(t_scale, v_scale, r_scale, t0);
            multipleIntegrate(y0, new_sol, arcs, 10);

            //for(int k = 0; k < y0.Length; k++)
            //    Debug.Log("y0[" + k + "] = " + y0[k]);

            if (insertedCoast)
            {
                if ( new_sol.tgo(new_sol.t0, arcs.Count-2) < 1 )
                {
                    // coast is less than one second, delete it and reconverge
                    RemoveArc(arcs, arcs.Count-2, new_sol);

                    if ( !runOptimizer(arcs) )
                    {
                        Fatal("failed to converge after removing zero length coast");
                        return;
                    }

                    new_sol = new Solution(t_scale, v_scale, r_scale, t0);
                    multipleIntegrate(y0, new_sol, arcs, 10);
                }
            }

            //for(int k = 0; k < y0.Length; k++)
                //Debug.Log("new y0[" + k + "] = " + y0[k]);

            this.solution = new_sol;
            //Debug.Log("done with bootstrap");

            yf = new double[arcs.Count*13];
            multipleIntegrate(y0, yf, arcs);

            //for(int k = 0; k < yf.Length; k++)
                //Debug.Log("new yf[" + k + "] = " + yf[k]);

        }

        public override void Optimize(double t0)
        {
            base.Optimize(t0);
            //Debug.Log("r_scale = " + r_scale);
            //Debug.Log("v_scale = " + v_scale);
            //Debug.Log("rTm = " + rTm * r_scale + " " + rTm);
            //Debug.Log("vTm = " + vTm * v_scale + " " + vTm);
        }
    }
}
