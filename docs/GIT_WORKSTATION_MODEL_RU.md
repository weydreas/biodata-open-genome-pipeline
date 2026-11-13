# Git-native модель рабочих станций Dante 9.1.0

## Каноническая версия

`main` в `weydreas/biodata-open-genome-pipeline` является единственной финальной production-веткой. Домашний компьютер, ноутбук и любые последующие рабочие станции являются равноправными Git checkout одного и того же репозитория.

## Поведение `start.cmd`

### Если активна `main`

Bootstrap:

1. определяет, что корень pipeline уже является Git checkout;
2. выполняет `git fetch origin main`;
3. если working tree чистый и локальный `main` является предком `origin/main`, выполняет только `git merge --ff-only origin/main`;
4. никогда не выполняет destructive `reset --hard` в Git-workstation checkout;
5. если есть локальные изменения, unpublished commits или divergence, ничего не переписывает и продолжает с локальным кодом с явным предупреждением.

Таким образом, после merge новой версии в GitHub следующий обычный запуск `start.cmd` на любой чистой станции автоматически подхватывает текущий `main`.

### Если активна feature/fix/docs ветка

Bootstrap распознаёт development checkout и запускает код текущей ветки без развёртывания `origin/main` поверх неё. Никакие environment-переменные для защиты рабочей ветки не требуются.

## Frozen launchers

`start.cmd` и `bootstrap_git.cmd` считаются frozen launchers. Git-native auto-fast-forward блокируется, если между локальным `HEAD` и `origin/main` изменился любой из этих файлов. Это исключает замену batch-файла во время его исполнения. Редкая миграция самого launcher выполняется отдельно после закрытия launcher-процесса.

## `.pipeline_repo`

`.pipeline_repo` больше не используется в нормальной Git-workstation установке. Он сохраняется только как backward-compatible механизм обновления старых runtime-инсталляций, у которых корень pipeline не является Git repository.

## Reference bundles

Startup больше не устанавливает open/research reference bundles автоматически. Он только проверяет их READY-статус. Установка или обновление выполняются явно через соответствующие пункты меню (`22` и `23`) тогда, когда эти данные действительно нужны для выбранной аналитической задачи.

## Локальные patient/runtime данные

`patients/*`, reference/runtime artifacts, `.active_patient.json`, secrets и тяжёлые геномные данные не входят в Git. Выбор клиента не должен делать working tree грязным и не должен блокировать безопасный auto-fast-forward `main`.

## Regression contract

Windows CI исполняет `tests/Test-GitNativeBootstrap.ps1`, который на синтетических локальных репозиториях доказывает:

- automatic fast-forward чистого `main`;
- отсутствие `.pipeline_repo` в Git checkout;
- блокировку live-update frozen launcher;
- восстановление auto-update после устранения launcher delta;
- сохранение активной feature-ветки без подмешивания `origin/main`.
