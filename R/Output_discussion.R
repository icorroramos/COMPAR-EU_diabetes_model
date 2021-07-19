### CREATE OUTPUTS FOR INTERNAL DISCUSSION WITH REGARD TO OUTPUTS IN PLATFORM

# Needs model results as created by 'Analysis script template.R' in global environment


# Create deltas for each cycle
f.int <- sim.results.female[[1]] %>%
  group_by(SDURATION) %>% 
  summarize_all(mean)

f.comp <- sim.results.female.comp[[1]] %>% 
  group_by(SDURATION) %>% 
  summarize_all(mean)

m.int <- sim.results.male[[1]] %>%
  group_by(SDURATION) %>% 
  summarize_all(mean)

m.comp <- sim.results.male.comp[[1]] %>% 
  group_by(SDURATION) %>% 
  summarize_all(mean)


f.int.ce <- data.frame(
  year = 1:nrow(f.int),
  Sex = 'Female',
  Treatment = 'SMI', 
  QALY = f.int$QALY,
  QALY_cum = cumsum(f.int$QALY),
  Med.cost = rowSums(f.int[ , 51:59]),
  Med.cost.cum = cumsum(rowSums(f.int[ , 51:59])),
  Tot.cost = f.int$TOTAL.COST
)

f.comp.ce <- data.frame(
  year = 1:nrow(f.comp),
  Sex = 'Female',
  Treatment = 'Usual care', 
  QALY = f.comp$QALY,
  QALY_cum = cumsum(f.comp$QALY),
  Med.cost = rowSums(f.comp[ , 51:59]),
  Med.cost.cum = cumsum(rowSums(f.comp[ , 51:59])),
  Tot.cost = f.comp$TOTAL.COST
)

m.int.ce<- data.frame(
  year = 1:nrow(m.int),
  Sex = 'Male',
  Treatment = 'SMI', 
  QALY = m.int$QALY,
  QALY_cum = cumsum(m.int$QALY),
  Med.cost = rowSums(m.int[ , 51:59]),
  Med.cost.cum = cumsum(rowSums(m.int[ , 51:59])),
  Tot.cost = m.int$TOTAL.COST
)

m.comp.ce <- data.frame(
  year = 1:nrow(m.comp),
  Sex = 'Male',
  Treatment = 'Usual care', 
  QALY = m.comp$QALY,
  QALY_cum = cumsum(m.comp$QALY),
  Med.cost = rowSums(m.comp[ , 51:59]),
  Med.cost.cum = cumsum(rowSums(m.comp[ , 51:59])),
  Tot.cost = m.comp$TOTAL.COST
)

ce.ot <- rbind(f.int.ce, f.comp.ce, m.int.ce, m.comp.ce)

# Plot for QALY
ggplot(ce.ot, aes(x = year, y = QALY, group = interaction(Sex, Treatment), color = Sex, linetype = Treatment)) +
  geom_line(size = 1) + 
  labs(x = 'Year since start of SMI', y = 'QALY gained per year')
  
ggplot(ce.ot, aes(x = year, y = QALY_cum, group = interaction(Sex, Treatment), color = Sex, linetype = Treatment)) +
  geom_line(size = 1) + 
  labs(x = 'Year since start of SMI', y = 'Total QALYs accumulated since start of SMI')

# Plot for costs
ggplot(ce.ot, aes(x = year, y = Med.cost, group = interaction(Sex, Treatment), color = Sex, linetype = Treatment)) +
  geom_line(size = 1) + 
  labs(x = 'Year since start of SMI', y = 'Medical costs per year')

ggplot(ce.ot, aes(x = year, y = Med.cost.cum, group = interaction(Sex, Treatment), color = Sex, linetype = Treatment)) +
  geom_line(size = 1) + 
  labs(x = 'Year since start of SMI', y = 'Total medical costs since start of SMI')


# Total cost summary data
cost.output <- rbind(sim_CE_results_female_table_comp, sim_CE_results_female_table, sim_CE_results_male_table_comp, sim_CE_results_male_table)
cost.ids <- data.frame(
  Sex = c('Female', 'Female', 'Male', 'Male'),
  Treatment = rep(c('Comparator', 'Intervention'), 2)
) 
cost.data.wide <- as.data.frame(cbind(cost.ids, cost.output[, 1:7]))

cost.data <- cost.data.wide %>% 
  gather(cost.type, value, 3:9)

med.cost.labs <- colnames(cost.data.wide[, 3:6])

med.cost.data <- cost.data %>% 
  filter(cost.type %in% med.cost.labs)

ggplot(cost.data, aes(fill = cost.type, y = value, x = Treatment)) +
  geom_bar(position = 'stack', stat = 'identity')

ggplot(med.cost.data, aes(fill = cost.type, y = value, x = Treatment)) +
  geom_bar(position = 'stack', stat = 'identity')
