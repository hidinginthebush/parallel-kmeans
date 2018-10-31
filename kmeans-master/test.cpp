#pragma once

#include <stdio.h>
#include <tchar.h>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <omp.h>

using namespace std;

struct Point
{
	int ClusterID;
	double CoOrdinat1;
	double CoOrdinat2;
	double CoOrdinat3;
	double CoOrdinat4;

	double operator-(Point const& rhs) {
		double tmp = 0;

		tmp += abs(rhs.CoOrdinat1 - this->CoOrdinat1) +
			abs(rhs.CoOrdinat2 - this->CoOrdinat2) +
			abs(rhs.CoOrdinat3 - this->CoOrdinat3) +
			abs(rhs.CoOrdinat4 - this->CoOrdinat4);
		return tmp;
	}

public:
	Point() {}
	Point(int ClusterID, double CoOrdinat1, double CoOrdinat2, double CoOrdinat3, double CoOrdinat4)
	{
		this->ClusterID = ClusterID;
		this->CoOrdinat1 = CoOrdinat1;
		this->CoOrdinat2 = CoOrdinat2;
		this->CoOrdinat3 = CoOrdinat3;
		this->CoOrdinat4 = CoOrdinat4;
	}
};

int NUM_CLUSTERS;
int NUM_POINTS;

//=======================Sequential Implementation==========================
vector<string> seperate(string line, char sep) {
	vector<string> data;
	int ind = 0;
	string tmp = "";

	for (int i = 0; i < line.length(); i++)
	{
		if (line[i] == sep) {
			data.push_back(tmp);
			tmp = "";
		}
		else
		{
			tmp += line[i];
		}
	}
	data.push_back(tmp);
	return data;
}
vector<Point> readData()
{
	ifstream IrisData;
	string data;
	vector<Point>points;
	// Opening Iris file To read data from
	IrisData.open("IrisDataset.txt");
	int c = 0;
	// Getting seperated data
	while (IrisData >> data) {
		Point myPoint;
		// Reading Points
		if (c == 0)
		{
			NUM_POINTS = stoi(data);
			c++;
			continue;
		}
		// Reading Clusters
		if (c == 1)
		{
			NUM_CLUSTERS = stoi(data);
			c++;
			continue;
		}

		vector<string> sepData = seperate(data, ',');
		myPoint.CoOrdinat1 = stod(sepData[0]);
		myPoint.CoOrdinat2 = stod(sepData[1]);
		myPoint.CoOrdinat3 = stod(sepData[2]);
		myPoint.CoOrdinat4 = stod(sepData[3]);
		points.push_back(myPoint);
		c++;
	}

	// closing file
	IrisData.close();

	return points;
};
vector<Point> cluster(vector<Point> clusters, vector<Point> data)
{
	for (int i = 0; i < NUM_POINTS; i++)
	{
		int C_ID = -1;
		double Min = 10000.0;
		for (int j = 0; j < NUM_CLUSTERS; j++)
		{
			double distance = sqrtl(powl((data[i].CoOrdinat1 - clusters[j].CoOrdinat1), 2) +
				powl((data[i].CoOrdinat2 - clusters[j].CoOrdinat2), 2) +
				powl((data[i].CoOrdinat3 - clusters[j].CoOrdinat3), 2) +
				powl((data[i].CoOrdinat4 - clusters[j].CoOrdinat4), 2));
			if (Min > distance)
			{
				Min = distance;
				C_ID = j;
			}
		}

		data[i].ClusterID = C_ID;
	}
	return data;
}
vector<Point> get_new_clusters(vector<Point> data) {

	vector<vector<Point>> data_categorized(NUM_CLUSTERS);

	for (int i = 0; i < data.size(); i++)
	{
		data_categorized[data[i].ClusterID].push_back(data[i]);
	}

	vector<Point> ReallocatedClusters(NUM_CLUSTERS);

	for (int i = 0; i < NUM_CLUSTERS; i++)
	{
		double tmp1 = 0, tmp2 = 0, tmp3 = 0, tmp4 = 0;
		for (int j = 0; j < data_categorized[i].size(); j++)
		{
			tmp1 += data_categorized[i][j].CoOrdinat1;
			tmp2 += data_categorized[i][j].CoOrdinat2;
			tmp3 += data_categorized[i][j].CoOrdinat3;
			tmp4 += data_categorized[i][j].CoOrdinat4;
		}

		Point tmp(i, (tmp1 / data_categorized[i].size()), (tmp2 / data_categorized[i].size()), (tmp3 / data_categorized[i].size()), (tmp4 / data_categorized[i].size()));
		ReallocatedClusters[i] = tmp;
	}
	return ReallocatedClusters;
}
bool convergence(vector<Point> old_centroids, vector<Point> clusters) {
	bool threshold = true;

	for (int i = 0; i < NUM_CLUSTERS; i++)
	{
		if (old_centroids[i] - clusters[i] > 0.001)
			threshold = false;
	}

	if (threshold)
		return true;

	for (int i = 0; i < NUM_CLUSTERS; i++)
	{
		if (old_centroids[i].CoOrdinat1 != clusters[i].CoOrdinat1 ||
			old_centroids[i].CoOrdinat2 != clusters[i].CoOrdinat2 ||
			old_centroids[i].CoOrdinat3 != clusters[i].CoOrdinat3 ||
			old_centroids[i].CoOrdinat4 != clusters[i].CoOrdinat4)
		{
			return false;
		}
	}

	return true;
}

//=======================Parallelized Implementation========================
vector<Point> cluster_P(vector<Point> clusters, vector<Point> data)
{
	#pragma omp parallel for schedule(guided)
		for (int i = 0; i < NUM_POINTS; i++)
		{
			int C_ID = -1;
			double Min = 10000.0;
			for (int j = 0; j < NUM_CLUSTERS; j++)
			{
				//printf("i: %d    num: %d\n", i, omp_get_thread_num());
				double distance = sqrtl(powl((data[i].CoOrdinat1 - clusters[j].CoOrdinat1), 2) +
					powl((data[i].CoOrdinat2 - clusters[j].CoOrdinat2), 2) +
					powl((data[i].CoOrdinat3 - clusters[j].CoOrdinat3), 2) +
					powl((data[i].CoOrdinat4 - clusters[j].CoOrdinat4), 2));
				if (Min > distance)
				{
					Min = distance;
					C_ID = j;
				}
			}

			data[i].ClusterID = C_ID;
		}
	
	return data;
}
vector<Point> get_new_clusters_P(vector<Point> data) {
	int t = omp_get_num_threads();
	vector<vector<Point>> data_categorized(NUM_CLUSTERS);
	vector<Point> ReallocatedClusters(NUM_CLUSTERS);
	#pragma omp parallel shared(data_categorized, ReallocatedClusters)
	{
		#pragma omp for schedule(guided)
		for (int i = 0; i < data.size(); i++)
		{
			//t = omp_get_num_threads();
			#pragma omp critical 
			data_categorized[data[i].ClusterID].push_back(data[i]);
		}

		#pragma omp barrier 
		#pragma omp for schedule(guided)
		for (int i = 0; i < NUM_CLUSTERS; i++)
		{
			double tmp1 = 0, tmp2 = 0, tmp3 = 0, tmp4 = 0;
			for (int j = 0; j < data_categorized[i].size(); j++)
			{
				tmp1 += data_categorized[i][j].CoOrdinat1;
				tmp2 += data_categorized[i][j].CoOrdinat2;
				tmp3 += data_categorized[i][j].CoOrdinat3;
				tmp4 += data_categorized[i][j].CoOrdinat4;
			}

			Point tmp(i, (tmp1 / data_categorized[i].size()), (tmp2 / data_categorized[i].size()), (tmp3 / data_categorized[i].size()), (tmp4 / data_categorized[i].size()));
			ReallocatedClusters[i] = tmp;
		}
	}
	
	return ReallocatedClusters;
}
bool convergence_P(vector<Point> old_centroids, vector<Point> clusters) {
	bool threshold = true;

#pragma omp parallel for schedule (guided)
	for (int i = 0; i < NUM_CLUSTERS; i++)
	{
		if (old_centroids[i] - clusters[i] > 0.001)
			threshold = false;
	}

	if (threshold)
		return true;

	for (int i = 0; i < NUM_CLUSTERS; i++)
	{
		if (old_centroids[i].CoOrdinat1 != clusters[i].CoOrdinat1 ||
			old_centroids[i].CoOrdinat2 != clusters[i].CoOrdinat2 ||
			old_centroids[i].CoOrdinat3 != clusters[i].CoOrdinat3 ||
			old_centroids[i].CoOrdinat4 != clusters[i].CoOrdinat4)
		{
			return false;
		}
	}

	return true;
	//} 



}


int main()
{

	#pragma region Sequential K-Means Algorithm
	
		// Reading Data from Iris file
		vector<Point> data_S = readData();
	
		// Intiating Clusters with arbitary points
		vector<Point> clusters_S(NUM_CLUSTERS);
		for (int i = 0; i < NUM_CLUSTERS; i++)
		{
			clusters_S[i] = data_S[i];
			clusters_S[i].ClusterID = i;
		}
	
		// Main Algorithm
		while (true)
		{
			vector<Point> old_centroids_S = clusters_S;
			// Assigning every point to a cluster
			data_S = cluster(clusters_S, data_S);
	
			// Calculating new centroids
			clusters_S = get_new_clusters(data_S);
	
			if (convergence(old_centroids_S, clusters_S))
				break;
		}
	
	
		// Write out Centroids to the output file
		ofstream file;
		file.open("dataset2.txt", ios_base::trunc);
		file << "\t\t This is the Iris data set Output file for " << NUM_CLUSTERS << " clusters\n\n";
		file << "Cluster #\tCoordinate 1\tCoordinate 2\tCoordinate 3\tCoordinate 4\n";
		for (int i = 0; i < NUM_CLUSTERS; i++)
		{
			file << clusters_S[i].ClusterID << " \t "
				<< clusters_S[i].CoOrdinat1 << " \t "
				<< clusters_S[i].CoOrdinat2 << " \t "
				<< clusters_S[i].CoOrdinat3 << " \t "
				<< clusters_S[i].CoOrdinat4 << endl;
		}
		file.close();
	#pragma endregion

	#pragma region Parallelized K-Means Algorithm

	omp_set_num_threads(20);
	// Reading Data from Iris file
	vector<Point> data = readData();

	// Intiating Clusters with arbitary points
	vector<Point> clusters(NUM_CLUSTERS);
	#pragma omp parallel 
	{
		#pragma omp for schedule(guided)
		for (int i = 0; i < NUM_CLUSTERS; i++)
		{
			clusters[i] = data[i];
			clusters[i].ClusterID = i;
			//printf("i: %d    num: %d\n", i, omp_get_thread_num());
		}
	}
	

	// Main Algorithm
	while (true)
	{
		vector<Point> old_centroids = clusters;
		// Assigning every point to a cluster
		data = cluster_P(clusters, data);

		// Calculating new centroids
		clusters = get_new_clusters_P(data);

		if (convergence_P(old_centroids, clusters))
			break;
	}


	// Write out Centroids to the output file
	file.open("dataset2.txt", ios_base::trunc);
	file << "\t\t This is the Iris data set Output file for " << NUM_CLUSTERS << " clusters\n\n";
	file << "Cluster #\tCoordinate 1\tCoordinate 2\tCoordinate 3\tCoordinate 4\n";

	#pragma omp parallel for schedule(guided)
	for (int i = 0; i < NUM_CLUSTERS; i++)
	{
		#pragma omp critical
		file << clusters[i].ClusterID << " \t "
			<< clusters[i].CoOrdinat1 << " \t "
			<< clusters[i].CoOrdinat2 << " \t "
			<< clusters[i].CoOrdinat3 << " \t "
			<< clusters[i].CoOrdinat4 << endl;
	}

	file.close();
#pragma endregion


	return 0;
}
