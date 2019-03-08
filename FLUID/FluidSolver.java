

/**
 * Fluid solver with vorticity confinement and buoyancy force.
 *
 * Used for teaching purpose at ISTY 
 *
 * Based on original work from Alexander McKenzie, Caltech
 * This Class is basically sending in/out data from the Webserver to the Core computing C libraries
 **/

public class FluidSolver
{
    int n, size;
    float dt;
    float visc = 0.0f;
    float diff = 0.0f;
    float[] d, dOld;
    float[] u, uOld;
    float[] v, vOld;
    float[] curl;

    /**
     * Set the grid size and timestep.
     **/
    public void setup(int n, float dt)
    {
        this.n = n;
        this.dt = dt;
        size = (n + 2) * (n + 2);

        reset();
    }

    /**
     * Reset the datastructures.
     * We use 1d arrays for speed.
     **/
    public void reset()
    {
        d    = new float[size];
        dOld = new float[size];
        u    = new float[size];
        uOld = new float[size];
        v    = new float[size];
        vOld = new float[size];
        curl = new float[size];

        for (int i = 0; i < size; i++)
        {
            u[i] = uOld[i] = v[i] = vOld[i] = 0.0f;
            d[i] = dOld[i] = curl[i] = 0.0f;
        }

    }

    /**
     * The basic velocity solving routine as described by Stam. This function is simply a stub to call the C routine
     **/
    public void velocitySolver()
    {
    fluid.c_velocitySolver(u, uOld, v, vOld, curl, d, visc, dt, n, size);
    }

    /**
     * The basic density solving routine. This function is simply a stub to call the C routine
     **/
    public void densitySolver()
    {
    	//System.out.println("visc=" + Arrays.toString(dOld)); 

    // add density inputted by mouse
    fluid.c_densitySolver(d, dOld, diff, u, v , dt, n, size);
    }
}
