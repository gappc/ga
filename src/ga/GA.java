package ga;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class GA {

	private int POPULATION_SIZE = 50;
	private int GENOME_LENGTH = 5;
	
	private int ITERATIONS = (int)1e6;

	private int MAX_INIT_VALUE = 600;
	private int MIN_INIT_VALUE = -600;
	
	private double[] values;
	
	public static void main(String[] args) {
		GA ga = new GA();
		ga.iteration();
	}

	public double[] iteration() {
		Random rand = new Random();

		boolean end = false;
		int counter = 0;

		// Generate random population
		List<double[]> population = initPopulation();

		while (!end) {
			// recombination
			List<double[]> offsprings = new ArrayList<double[]>();
			for (int i = 0; i < POPULATION_SIZE; i++) {
				int indexParent1 = rand.nextInt(POPULATION_SIZE);
				int indexParent2 = rand.nextInt(POPULATION_SIZE);

				offsprings.add(recombination(population.get(indexParent1),
						population.get(indexParent2)));
			}
			
			// mutation
			List<double[]> mutated = new ArrayList<double[]>();
			for (int i = 0; i < POPULATION_SIZE; i++) {
				mutated.add(mutation(offsprings.get(i)));
			}
			
			// evaluation
			values = new double[POPULATION_SIZE * 2];
			for (int i = 0; i < POPULATION_SIZE; i++) {
				values[i] = fitness(population.get(i));
			}
			for (int i = 0; i < POPULATION_SIZE; i++) {
				values[i + POPULATION_SIZE] = fitness(mutated.get(i));
			}
			
			// selection
			for (int i = 0; i < POPULATION_SIZE * 2; i++) {
				for (int j = i + 1; j < POPULATION_SIZE * 2; j++) {
					if (values[j] < values[i]) {
						double tmp = values[j];
						values[j] = values[i];
						values[i] = tmp;
						
						double[] tmpGenome;
						if (i < POPULATION_SIZE && j < POPULATION_SIZE) {
							//i and j refer to population
							tmpGenome = population.get(i);
							population.set(i, population.get(j));
							population.set(j, tmpGenome);
						}
						else if (i >= POPULATION_SIZE && j < POPULATION_SIZE) { 
							//i refer to mutated, j refer to population
							tmpGenome = mutated.get(i - POPULATION_SIZE);
							mutated.set(i - POPULATION_SIZE, population.get(j));
							population.set(j, tmpGenome);
						}
						else if (i < POPULATION_SIZE && j >= POPULATION_SIZE) { 
							//i refer to population, j refer to mutated
							tmpGenome = population.get(i);
							population.set(i, mutated.get(j - POPULATION_SIZE));
							mutated.set(j - POPULATION_SIZE, tmpGenome);
						}
						else {
							//i and j refer to mutated
							tmpGenome = mutated.get(i - POPULATION_SIZE);
							mutated.set(i - POPULATION_SIZE, mutated.get(j - POPULATION_SIZE));
							mutated.set(j - POPULATION_SIZE, tmpGenome);
						}
					}
				}
			}
			counter++;
			if (counter == ITERATIONS) {
				end = true;
			}
			
			if (counter % 1e4 == 0 || counter < 10) {
				System.out.println(counter + " " + values[0]);
				printGenome(population.get(0));
			}
		}
		System.out.println("Fitness value: " + values[0]);
		printGenome(population.get(0));
		return null;
	}

	private List<double[]> initPopulation() {
		List<double[]> population = new ArrayList<double[]>();
		Random rand = new Random();
		int intervall = MAX_INIT_VALUE - MIN_INIT_VALUE;
		
		for (int i = 0; i < POPULATION_SIZE; i++) {
			double[] singlePopulation = new double[GENOME_LENGTH];
			for (int j = 0; j < GENOME_LENGTH; j++) {
				singlePopulation[j] = rand.nextDouble() * intervall - intervall / 2;
			}
			population.add(singlePopulation);
		}
		return population;
	}

	private double[] recombination(double[] parent1, double[] parent2) {
		assert (parent1.length == parent2.length);
		double[] recombined = new double[parent1.length];
		for (int i = 0; i < parent1.length; i++) {
			recombined[i] = (parent1[i] + parent2[i]) / 2;
		}
		return recombined;
	}

	private double[] mutation(double[] genome) {
		double[] mutatedPopulation = new double[genome.length];

		double mutationProbability = 1.0 / genome.length;
		Random rand = new Random();

		for (int i = 0; i < genome.length; i++) {
			mutatedPopulation[i] = genome[i];
			if (rand.nextDouble() < mutationProbability) {
				if (rand.nextDouble() < 0.5) {
					mutatedPopulation[i] += rand.nextDouble();
				} else {
					mutatedPopulation[i] -= rand.nextDouble();
				}
			}
		}
		return mutatedPopulation;
	}
	
	private double fitness(double[] solution) {
		double result = 1;
		double sum = 0;
		double prod = 1;
		for (int i = 0; i < solution.length; i++) {
			sum += solution[i] * solution[i];
			prod *= Math.cos(solution[i] / Math.sqrt(i + 1));
		}

		result += sum / 4000 - prod;
		return result;
	}
	
	private void printGenome(double[] solution) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < GENOME_LENGTH; i++) {
			sb.append(solution[i] + " | ");
		}
		System.out.println(sb.toString());
	}
}
