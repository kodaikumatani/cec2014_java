import java.lang.Math;
import java.util.Random;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
/* csv */
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.IOException;

public class ES {
	public static int f = 1;	// 0, 1, 2, 3, ..., 22
	public static int dim = 10;	// 10 (or 30, 50, 100)
	public static double MIN = -100.0;
	public static double MAX = 100.0;
	
	public static void main(String[] args) {

		int ramda = 10;
		int mu = 2;
		int rand = 1;
		int trial = 11;
		
		//従来手法
		System.out.println(calculation(100000, ramda, mu, rand, false));
		experiment(false, trial);

		//改良手法１
		System.out.println(calculation(100000, ramda, mu, rand, true));
		experiment(true, trial);

		//改良手法２
		double E1 = calculation(50000, ramda, mu, rand, false);
		double E2 = calculation(50000, ramda, mu, rand, true);
		System.out.printf("Fixed:%f, Annealing:%f%n", E1, E2);
		experiment_imp2(trial);

		//改良手法３
		System.out.println(imp3Calculation(100000));
		experiment_imp3(trial);

		/*
		double dist = Landscape.distance(f, x);
		System.out.println("Distance = " + dist);
		double[] fval = new double[Landscape.numpoints+1];
		fval = Landscape.evaluate(f, x);
		for(int i=0; i<=Landscape.numpoints; i++){
			System.out.println(fval[i]);
		}
		*/
	}
	

	private static void experiment(boolean annealing, int trial)	{
		List<Double[]> data = new ArrayList<>();
		int eval_times = 100000;
		List<Integer> ramda = Arrays.asList(10, 20, 40, 80, 100);
		List<Double> mu = Arrays.asList(0.2, 0.4, 0.6, 0.8);
		List<Integer> rand = Arrays.asList(1, 5, 10, 15);

		int count = 1;
		for (Integer i : ramda) {
			for (Double j : mu) {
				for (Integer k : rand) {
					Double[] eval = new Double[trial];
					for (int l = 0; l < trial; l++) {
						eval[l] = calculation(eval_times, i, (int)(i*j), k, annealing);
					}
					Double[] cpy = Arrays.copyOf(eval,trial);
					Arrays.sort(cpy);
					double median = cpy[cpy.length / 2];
					
					Double[] line = new Double[trial+4];
					System.arraycopy(new Double[]{(double)i, j, (double)k, median}, 0, line, 0, 4);
					System.arraycopy(eval, 0, line, 4, trial);

					data.add(line);
					progressBar(count++, ramda.size() * mu.size() * rand.size());
				}
			}
		}
		try {
			FileWriter fw;
			if (annealing)	{
				fw = new FileWriter("data/imp1/func_"+Integer.valueOf(f).toString()+".csv", false);
			} else {
				fw = new FileWriter("data/conv/func_"+Integer.valueOf(f).toString()+".csv", false);
			}
            PrintWriter pw = new PrintWriter(new BufferedWriter(fw));
			
            for(Double[] line : data) {
				pw.print(line[0]+","+line[1]+","+line[2]+","+line[3]+",");
				for (int i = 0; i < line.length - 4; i++)	{
					pw.print(line[i+4]);
					pw.print(",");
				}
				pw.println();
            }
            pw.close();
 
        } catch (IOException ex) {
			System.out.println("Failed to save file");
            ex.printStackTrace();
        }
	}

	private static void experiment_imp2(int trial)	{
		List<Double[]> data = new ArrayList<>();
		int eval_times = 50000;
		List<Integer> ramda = Arrays.asList(10, 20, 40, 80, 100);
		List<Double> mu = Arrays.asList(0.2, 0.4, 0.6, 0.8);
		List<Integer> rand = Arrays.asList(1, 5, 10, 15);

		int count = 1;
		for (Integer i : ramda) {
			for (Double j : mu) {
				for (Integer k : rand) {
					Double[] eval = new Double[trial];
					Double[] fixed = new Double[trial];
					Double[] annealing = new Double[trial];
					for (int l = 0; l < trial; l++) {
						fixed[l] = calculation(eval_times, i, (int)(i*j), k, false);
						annealing[l] = calculation(eval_times, i, (int)(i*j), k, true);
						eval[l] = fixed[l] < annealing[l] ? fixed[l] : annealing[l];
					}
					Double[] cpy = Arrays.copyOf(eval,trial);
					Arrays.sort(cpy);
					double median = cpy[cpy.length / 2];
					
					Double[] line = new Double[trial*3+4];
					System.arraycopy(new Double[]{(double)i, j, (double)k, median}, 0, line, 0, 4);
					System.arraycopy(eval, 0, line, 4, trial);
					System.arraycopy(fixed, 0, line, 4+trial, trial);
					System.arraycopy(annealing, 0, line, 4+trial*2, trial);

					data.add(line);
					progressBar(count++, ramda.size() * mu.size() * rand.size());
				}
			}
		}
		try {
			FileWriter fw = new FileWriter("data/imp2/func_"+Integer.valueOf(f).toString()+".csv", false);
            PrintWriter pw = new PrintWriter(new BufferedWriter(fw));
			
            for(Double[] line : data) {
				pw.print(line[0]+","+line[1]+","+line[2]+","+line[3]+",");
				for (int i = 0; i < line.length - 4; i++)	{
					pw.print(line[i+4]);
					pw.print(",");
				}
				pw.println();
            }
            pw.close();
 
        } catch (IOException ex) {
			System.out.println("Failed to save file");
            ex.printStackTrace();
        }
	}

	private static void experiment_imp3(int trial)	{
		List<Double[]> data = new ArrayList<>();
		int eval_times = 100000;

		Double[] eval = new Double[trial];
		for (int i = 0; i < trial; i++) {
			eval[i] = imp3Calculation(100000);
		}
		Double[] median = Arrays.copyOf(eval,trial);
		Arrays.sort(median);
		try {
			FileWriter fw = new FileWriter("data/imp3.csv", true);
            PrintWriter pw = new PrintWriter(new BufferedWriter(fw));
			
			pw.print("F"+Integer.valueOf(f).toString()+",");
			pw.print(median[median.length / 2]+",");
			for (int i = 0; i < trial; i++)	{
				pw.print(eval[i]);
				pw.print(",");
			}
			pw.println();
            pw.close();
 
        } catch (IOException ex) {
			System.out.println("Failed to save file");
            ex.printStackTrace();
        }
	}

	private static double calculation(int eval_times, int ramda, int mu, double rand, boolean annealing)	{
		List<Double[]> data = new ArrayList<>();
		boolean out = false;
		int max_gene = (int) eval_times / ramda;

		Map<Double, double[]> children = getFirstChildren(ramda);
		Map<Double, double[]> parents = selection(mu, children);
		if (out) data.add(children.keySet().toArray(new Double[children.size()]));
		for (int i = 1; i < max_gene; i++)	{
			children = getChildren(i, max_gene, ramda, parents, rand, annealing);
			parents = selection(mu, concatenate(parents, children));
			if (out) data.add(children.keySet().toArray(new Double[children.size()]));
		}
		if (out) exportCsv(data);
		return parents.keySet().toArray(new Double[mu])[0];
	}

	private static double imp3Calculation(int eval_times)	{
		List<Double[]> data = new ArrayList<>();
		boolean out = false;
		double sync = 0.5;
		int ramda = 100;
		int mu = 80;
		double rand = 20;
		int first = (int)(eval_times * sync) / ramda;

		Map<Double, double[]> children = getFirstChildren(ramda);
		Map<Double, double[]> parents = selection(mu, children);
		if (out) data.add(children.keySet().toArray(new Double[children.size()]));
 		for (int i = 1; i < first; i++)	{
			children = getChildren(i, first, ramda, parents, rand, true);
			parents = selection(mu, concatenate(parents, children));
			if (out) data.add(children.keySet().toArray(new Double[children.size()]));
		}

		//スイッチ
		ramda = 5;
		mu = 1;
		rand = 3;
		parents = selection(mu, parents);
		
		int second = (int)(eval_times * (1-sync)) / ramda;
		for (int i = 0; i < second; i++)	{
			children = getChildren(i, second, ramda, parents, rand, true);
			parents = selection(mu, concatenate(parents, children));
			if (out) data.add(children.keySet().toArray(new Double[children.size()]));
		}
		if (out)  exportCsv(data);

		return parents.keySet().toArray(new Double[mu])[0];
	}

	public static Map<Double, double[]> getFirstChildren(int ramda) {
		Map<Double, double[]> children = new TreeMap<>();
		for (int i = 0; i < ramda; i++) {
			double[] vector = new Random().doubles(dim, MIN, MAX).toArray();
			double value = Simulator.evaluate(f, vector);
			children.put(value, vector);
		}
		return children;
	}

	public static Map<Double, double[]> selection(int mu, Map<Double, double[]> map) {
		Map<Double, double[]> selected = new TreeMap<>();
		Set<Double> keys = map.keySet();
		for (int i = 0; i < mu; i++) {
			Double key = keys.toArray(new Double[0])[i];
			//System.out.println(key + " => " + map.get(key));
			selected.put(key,map.get(key));
		}
		return selected;
	}

	public static Map<Double, double[]> getChildren(int gene, int max_gene, int ramda, Map<Double, double[]> parents, double rand, boolean annealing) {		
		Map<Double, double[]> children = new TreeMap<>();
		for (int i = 0; i < ramda; i++)	{
			double[] step_size = generateStepsize(gene, max_gene, annealing, rand);
			// ランダムな親から子供を生成する
			Random generator = new Random();
			Double[] keys = parents.keySet().toArray(new Double[parents.size()]);
			double[] p_vector = parents.get(keys[generator.nextInt(keys.length)]);

			double[] vector = new double[dim];
			for(int j = 0; j < dim; j++)	{
				double coordinate = (double)(p_vector[j] + step_size[j]);
				coordinate = coordinate > MAX ? MAX : coordinate;
				coordinate = coordinate < MIN ? MIN : coordinate;
				vector[j] = coordinate;
			}
			double value = Simulator.evaluate(f, vector);
			children.put(value, vector);
		}
		//System.out.println(randomValue.getClass().getSimpleName());
		return children;
	}

	public static double[] generateStepsize(int gene, int max_gene, boolean annealing, double rand)	{
		double step_range;
		if (annealing) {
			step_range = rand - ((double)rand / max_gene) * gene;
		} else {
			step_range = rand;
		}
		return new Random().doubles(dim, -step_range, step_range).toArray();
	}

	private static Map<Double, double[]> concatenate(Map<Double, double[]> map_A, Map<Double, double[]> map_B) {
		Map<Double, double[]> map = new TreeMap<>();
		map_A.forEach((key, value) -> {
			map.put(key, value);
		});
		map_B.forEach((key, value) -> {
			map.put(key, value);
		});
		return map;
	}

	public static void exportCsv(List<Double[]> data)	{
        try {
			FileWriter fw = new FileWriter("scatter.csv", false);
			PrintWriter pw = new PrintWriter(new BufferedWriter(fw));
			for(Double[] line : data) {
				for (int i = 0; i < line.length; i++)	{
					pw.println(line[i]);
				}
			}
			pw.close();
		} catch (IOException ex) {
			System.out.println("Failed to save file");
			ex.printStackTrace();
		}
    }

	public static void progressBar(int i, int trial) {
        double d = (double)i/trial;
		String p = String.format("%6.2f", d * 100);
		double p_cnt = (d * 100) / 10;
		p_cnt *= 2; // 100%で20個表示
		StringBuffer buf = new StringBuffer();
		buf.append(p);
		buf.append("% [");
		for (int j = 0; j < 20; j++) {
			if (j < (int)p_cnt) {
				buf.append("=");
			} else {
				buf.append(" ");
			}
		}
		buf.append("]");
		buf.append(Integer.toString(i));
		if (i == trial) {
			buf.append("\n");
		} else {
			buf.append("\r");
		}
		System.out.print(buf.toString());
	}
}