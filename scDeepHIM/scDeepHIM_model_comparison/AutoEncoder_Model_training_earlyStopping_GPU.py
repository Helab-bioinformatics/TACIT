import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import os
import time
import matplotlib.pyplot as plt
from torch.cuda import memory_allocated, max_memory_allocated, reset_max_memory_allocated

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

class AutoEncoder(nn.Module):
    def __init__(self):
        super(AutoEncoder, self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(3 * 273119, 512),
            nn.BatchNorm1d(512),
            nn.LeakyReLU(0.1),
            nn.Dropout(0.1),
            nn.Linear(512, 256),
            nn.BatchNorm1d(256),
            nn.LeakyReLU(0.1),
            nn.Dropout(0.1),
            nn.Linear(256, 128),
            nn.BatchNorm1d(128),
            nn.LeakyReLU(0.1),
            nn.Dropout(0.1),
            nn.Linear(128, 64)
        )
        self.decoder = nn.Sequential(
            nn.Linear(64, 128),
            nn.BatchNorm1d(128),
            nn.LeakyReLU(0.1),
            nn.Dropout(0.1),
            nn.Linear(128, 256),
            nn.BatchNorm1d(256),
            nn.LeakyReLU(0.1),
            nn.Dropout(0.1),
            nn.Linear(256, 512),
            nn.BatchNorm1d(512),
            nn.LeakyReLU(0.1),
            nn.Dropout(0.1),
            nn.Linear(512, 3 * 273119),
            nn.BatchNorm1d(3 * 273119),
            nn.LeakyReLU(0.1)
        )

    def forward(self, x):
        x = x.view(-1, 3 * 273119)
        x = self.encoder(x)
        x = self.decoder(x)
        x = x.view(-1, 3, 273119)
        return x


def load_data(data_paths, start_idx, end_idx):
    input_data, target_data = [], []
    for i in range(start_idx, end_idx):
        input_file = os.path.join(data_paths[0], f'input_matrix_{i}.txt')
        target_file = os.path.join(data_paths[1], f'output_matrix_{i}.txt')
        input_array = np.loadtxt(input_file)
        target_array = np.loadtxt(target_file)
        input_data.append(input_array)
        target_data.append(target_array)
    input_data = torch.tensor(np.array(input_data), dtype=torch.float32).to(device)
    target_data = torch.tensor(np.array(target_data), dtype=torch.float32).to(device)
    return input_data, target_data

data_paths = [
    '/path1/input_training',
    '/path2/output_training',
    '/path3/input_testing',
    '/path4/output_testing',
    '/path5/input_validation',
    '/path6/output_validation'
]

input_data, target_data = load_data(data_paths, 1, 128)
input_data1, target_data1 = load_data(data_paths[2:], 1, 32)
input_data2, target_data2 = load_data(data_paths[4:], 1, 6)


class EarlyStopping:
    def __init__(self, patience=50, delta=0.01):
        self.patience = patience
        self.delta = delta
        self.best_loss = np.inf
        self.counter = 0
        self.early_stop = False

    def __call__(self, loss):
        if self.best_loss - loss > self.delta:
            self.best_loss = loss
            self.counter = 0
        else:
            self.counter += 1
            if self.counter >= self.patience:
                self.early_stop = True


model = AutoEncoder().to(device)
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.005)


num_epochs = 1000
losses = []
correlations = []
test_losses = []
test_correlations = []
test2_losses = []
test2_correlations = []


gpu_memory_log = []
start_time = time.time()


early_stopping = EarlyStopping(patience=50, delta=0.01)


for epoch in range(num_epochs):
    model.train()
    outputs = model(input_data)
    loss = criterion(outputs, target_data)
    
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

    losses.append(loss.item())

    with torch.no_grad():
        output_array = outputs.cpu().numpy().flatten()
        target_array = target_data.cpu().numpy().flatten()
        correlation = np.corrcoef(output_array, target_array)[0, 1]
        correlations.append(correlation)


    with torch.no_grad():
        model.eval()
        test_outputs = model(input_data1)
        test_loss = criterion(test_outputs, target_data1)
        test_output_array = test_outputs.cpu().numpy().flatten()
        test_target_array = target_data1.cpu().numpy().flatten()
        test_correlation = np.corrcoef(test_output_array, test_target_array)[0, 1]
        test_losses.append(test_loss.item())
        test_correlations.append(test_correlation)


    with torch.no_grad():
        test_outputs2 = model(input_data2)
        test2_loss = criterion(test_outputs2, target_data2)
        test_output_array2 = test_outputs2.cpu().numpy().flatten()
        test_target_array2 = target_data2.cpu().numpy().flatten()
        test2_correlation = np.corrcoef(test_output_array2, test_target_array2)[0, 1]
        test2_losses.append(test2_loss.item())
        test2_correlations.append(test2_correlation)

    if (epoch + 1) % 1 == 0:
        current_memory = max_memory_allocated(device)
        gpu_memory_log.append(current_memory)
        reset_max_memory_allocated(device)
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}, Correlation: {correlation:.4f}, '
              f'Test Loss: {test_loss.item():.4f}, Test Correlation: {test_correlation:.4f}, '
              f'Test2 Loss: {test2_loss.item():.4f}, Test2 Correlation: {test2_correlation:.4f}, '
              f'GPU Memory: {current_memory / 1024**2:.2f} MB')


    early_stopping(loss.item())
    if early_stopping.early_stop:
        print("Early stopping triggered.")
        break

torch.save(model.state_dict(), 'autoencoder_model.pt')

end_time = time.time()
total_time = end_time - start_time

with open('gpu_memory_time_log.txt', 'w') as f:
    f.write(f"Total Training Time: {total_time:.2f} seconds\n")
    f.write("GPU Memory Usage per Epoch (MB):\n")
    f.write('\n'.join([f"{mem / 1024**2:.2f}" for mem in gpu_memory_log]))

fig, axs = plt.subplots(2, 1, figsize=(8, 12))

axs[0].plot(losses, label='Training Loss')
axs[0].plot(test_losses, label='Testing Loss')
axs[0].plot(test2_losses, label='Testing embryos Loss')
axs[0].set_xlabel('Epoch')
axs[0].set_ylabel('Loss')
axs[0].set_title('Training Loss')
axs[0].legend()

axs[1].plot(correlations, label='Correlation')
axs[1].plot(test_correlations, label='Test Correlation')
axs[1].plot(test2_correlations, label='Testing embryos Correlations')
axs[1].set_xlabel('Epoch')
axs[1].set_ylabel('Correlation')
axs[1].set_title('Correlation')
axs[1].legend()

plt.tight_layout()
plt.savefig('autoencoder_training_plots.svg')
plt.show()
